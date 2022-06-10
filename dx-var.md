### Hardware and OS

Acquired a [Dell Precision 3640](https://www.dell.com/en-us/work/shop/desktops-all-in-one-pcs/precision-3640-tower-workstation/spd/precision-3640-workstation) tower workstation with the following specs:

- Intel Xeon W-1290P (best [single-thread performance](https://www.cpubenchmark.net/singleThread.html#server-thread) in early 2021)
- 128GB DDR4-2666 ECC Memory (64GB minimum; ECC reduces odds of data corruption)
- 2x 512GB Toshiba Kioxia XG6 PCIe NVMe SSD (nvme0n1, nvme1n1)
- 2x 8TB 7200rpm SATA 3.5" HDD (sda, sdb)

1. On a Windows PC, download the [Ubuntu Desktop 20.04.4 ISO](https://mirrors.ocf.berkeley.edu/ubuntu-releases/20.04.4) image.
2. Follow [these instructions](https://ubuntu.com/tutorials/create-a-usb-stick-on-windows) to create a bootable USB stick.
3. Boot into BIOS on the new workstation (press F2 during Dell logo), and ensure that:
    - `SATA Operation` is set to AHCI (disables hardware RAID)
    - `Secure Boot` is enabled (works with Ubuntu for now)
    - `AC Recovery` is set to `Last Power State`
    - `POST Behaviour` is set to `Continue on Warnings` and `Keyboard Error Detection` is disabled
    - `Virtualization`, `VT for Direct I/O`, and `Trusted Execution` are enabled
4. Plug in the USB stick and boot into it (press F12 during Dell logo).
5. Choose to `Install Ubuntu` and select `Minimal installation` when prompted.
6. Under `Installation Type` choose `Something else`, and partition the first SSD (nvme0n1) as follows:
    - Two partitions, 128MB for EFI and remainder `ext4` mounted as `/`
    - Make sure the boot loader is installed on this SSD (nvme0n1)
7. Follow on-screen instructions for remaining steps until reboot (hostname set to `dx-var`).
8. After booting into the installed OS, start Terminal and install the following:
    ```bash
    sudo apt update
    sudo apt install -y openssh-server zfsutils-linux samba git vim screen parallel tree curl
    ```
9. Create a unix group for your team (E.g. `dx` below) and make it your primary group:
    ```bash
    sudo groupadd dx
    sudo usermod -g dx $USER
    sudo groupdel $USER
    ```
10. Create a mirrored ZFS disk pool on the two HDDs, with an L2ARC cache on the second SSD, and make it writable for your team:
    ```bash
    sudo zpool create hot mirror /dev/sda /dev/sdb
    sudo chown $USER:dx /hot
    sudo chmod 775 /hot
    sudo zpool add hot cache nvme1n1
    ```
11. Share the folder `/hot` over the network via samba:
    ```bash
    echo -e "[hot]\n  path=/hot\n  writeable=yes" | sudo tee -a /etc/samba/smb.conf
    sudo service smbd restart
    ```
12. Set a samba password for your username using `sudo smbpasswd -a $USER`
13. Allow OpenSSH/Samba connections through the firewall and then enable it:
    ```bash
    sudo ufw allow OpenSSH
    sudo ufw allow Samba
    sudo ufw enable
    ```

### Run Output

Since we set up `/hot` as network-attached storage (NAS) on the server, we can create a folder `/hot/runs` and specify `\\dx-var\hot\runs` as the run output folder on each sequencer. This is sufficient for small clinical labs with at least one engineer who can be "on call" if anything breaks. However, Illumina strongly recommends an enterprise NAS solution hosted by your organization's IT team to ensure redundancy, reliability, and tech support. Our IT team configured a 10TB NAS (`\\aipnasvm01\pathnovaseq`) co-located with our Dragen server in a nearby data center, so we went with that. Below is how we mounted this NAS under `/mnt/pns` on our server.

1. Install dependencies, create the mount point, and a file to store NAS credentials:
    ```bash
    sudo apt install -y cifs-utils
    sudo mkdir /mnt/pns
    sudo chmod 775 /mnt/pns
    sudo touch /etc/cifs-utils/aipnasvm01-creds
    sudo chmod 600 /etc/cifs-utils/aipnasvm01-creds
    ```
2. Edit `/etc/cifs-utils/aipnasvm01-creds` and enter credentials in the following format:
    ```
    username=ckandoth
    password=blahblah
    domain=ad
    ```
3. Create `/etc/systemd/system/mnt-pns.mount` with content as follows:
    ```
    [Unit]
    Description=pathnovaseq share on aipnasvm01
    Requires=systemd-networkd.service
    After=network-online.target
    Wants=network-online.target

    [Mount]
    What=//aipnasvm01.ad.medctr.ucla.edu/pathnovaseq
    Where=/mnt/pns
    Options=credentials=/etc/cifs-utils/aipnasvm01-creds,iocharset=utf8,rw,x-systemd.automount,uid=1000,gid=1001,dir_mode=0775,file_mode=0664
    Type=cifs
    TimeoutSec=30

    [Install]
    WantedBy=multi-user.target
    ```
4. Start and enable the service, then run `df -ht cifs` to make sure the NAS is mounted:
    ```bash
    sudo systemctl start mnt-pns.mount
    sudo systemctl enable mnt-pns.mount
    ```
5. Create `/etc/systemd/system/mnt-pns.automount` with content as follows:
    ```
    [Unit]
    Description=pathnovaseq share on aipnasvm01

    [Automount]
    Where=/mnt/pns
    TimeoutIdleSec=0

    [Install]
    WantedBy=multi-user.target
    ```
6. Enable this service to ensure that the NAS is automounted on demand:
    ```bash
    sudo systemctl enable mnt-pns.automount
    ```

### Run Scanner

[Click here](https://miso-lims.readthedocs.io/projects/runscanner) to read about Run Scanner. These are the installation steps summarized:

1. Install dependencies, open a port, and create a config file:
    ```bash
    sudo apt install -y tomcat9 maven clang cmake libjsoncpp-dev autoconf libtool build-essential
    sudo ufw allow 8080/tcp
    echo -e '<Context>\n    <Parameter name="runscanner.configFile" value="/etc/runscanner.json" override="false"/>\n</Context>' | sudo tee /var/lib/tomcat9/conf/Catalina/localhost/runscanner.xml
    ```
2. Create `/etc/runscanner.json` and use it to describe paths to each run folder. For example:
    ```json
    [
        {
            "path": "/hot/runs",
            "platformType": "ILLUMINA",
            "name": "default",
            "timeZone": "America/Los_Angeles",
            "parameters": {}
        },
        {
            "path": "/mnt/pns/runs",
            "platformType": "ILLUMINA",
            "name": "default",
            "timeZone": "America/Los_Angeles",
            "parameters": {}
        }
    ]
    ```
3. Download, build, and deploy runscanner with support for parsing Illumina InterOp files:
    ```bash
    git clone --recursive git@github.com:miso-lims/runscanner.git
    cd runscanner/runscanner-illumina
    export CXX=$(which clang++) && ./build-illumina-interop && autoreconf -i && ./configure && make && sudo make install
    cd ..
    mvn clean package
    sudo cp scanner/target/scanner-*.war /var/lib/tomcat9/webapps/runscanner.war
    ```
4. Open `localhost:8080/runscanner` in a browser to make sure the graphical interface works.

### VarSeq

[Click here](https://www.goldenhelix.com/learning/article-categories/varseq-tutorials/) to see some tutorials on this suite of highly configurable tools. VSPipeline helps automate most of the steps from variant calling to a clinical report, and VSClinical provides a graphical interface for variant interpretation before sign-out by a molecular pathologist. These are the installation steps summarized:

1. Download the version of VarSeq we want to deploy, and install it:
    ```bash
    sudo mkdir /opt/VarSeq-2.2.5
    sudo chown $USER:dx /opt/VarSeq-2.2.5
    curl -L https://www.goldenhelix.com/download/VarSeq/VarSeq-Lin64-2.2.5.tar.gz -o /opt/VarSeq-2.2.5/VarSeq-Lin64-2.2.5.tar.gz
    tar -zxf /opt/VarSeq-2.2.5/VarSeq-Lin64-2.2.5.tar.gz -C /opt/VarSeq-2.2.5 --strip-components 1
    ```
2. Create `/hot/varseq` for app data and symlink to it from the installation folder:
    ```bash
    mkdir /hot/varseq
    chmod 775 /hot/varseq
    cd /opt/VarSeq-2.2.5
    ln -s /hot/varseq AppData
    ```
3. Start VarSeq
