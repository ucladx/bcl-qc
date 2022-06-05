### Hardware and OS

Acquired a [Dell Precision 3640](https://www.dell.com/en-us/work/shop/desktops-all-in-one-pcs/precision-3640-tower-workstation/spd/precision-3640-workstation) tower workstation with the following specs:

- Intel Xeon W-1290P (best [single-thread performance](https://www.cpubenchmark.net/singleThread.html#server-thread) in early 2021)
- 128GB DDR4-2666 ECC Memory (64GB minimum; ECC prevents genomic data corruption)
- 2x 512GB Toshiba Kioxia XG6 PCIe NVMe SSD (decent latency best bandwidth)
- 2x 8TB 7200rpm SATA 3.5" HDD (poor latency and bandwidth, but best capacity)

1. On a Windows PC, download the [Ubuntu Desktop 22.04 ISO](https://mirrors.ocf.berkeley.edu/ubuntu-releases/22.04) image.
2. Follow [these instructions](https://ubuntu.com/tutorials/create-a-usb-stick-on-windows) to create a bootable USB stick.
3. Boot into BIOS on the new workstation (press F2 during Dell logo), and ensure that:
    - `SATA Operation` is set to AHCI (disables hardware RAID)
    - `Secure Boot` is enabled (works with Ubuntu for now)
    - `AC Recovery` is set to `Last Power State`
    - `POST Behaviour` is set to `Continue on Warnings` and `Keyboard Error Detection` is disabled
    - `Virtualization`, `VT for Direct I/O`, and `Trusted Execution` are enabled
4. Plug in the USB stick and boot into it (press F12 during Dell logo).
5. Choose to `Install Ubuntu` and select `Minimal installation` when prompted.
6. Under `Installation Type` choose `Something else`, and partition the SSDs as follows:
    - 512GB SSD: Two partitions, 128MB for EFI and remainder `ext4` mounted as `/`
    - 512GB SSD: One partition, `ext4` mounted as `/srv`
    - Make sure the boot loader is installed on the first SSD (nvme0n1)
7. Follow on-screen instructions for remaining steps until reboot (hostname set to `dx-hpc`).
8. After booting into the installed OS, start Terminal and install the following:
    ```bash
    sudo apt update
    sudo apt install -y openssh-server zfsutils-linux samba git vim screen parallel tree
    ```
9. Create a unix group for your team (E.g. `dx` below) and make it your primary group:
    ```bash
    sudo groupadd dx
    sudo usermod -g dx $USER
    ```
10. Create a mirrored ZFS disk pool for redundant storage on the two 8TB HDDs, and make it writable for your team:
    ```bash
    sudo zpool create hot mirror /dev/sda /dev/sdb
    sudo chown $USER:dx /hot
    sudo chmod 775 /hot
    ```
11. Share the folder `/hot` over the network via samba:
    ```bash
    echo -e "[hot]\n  path=/hot\n  writeable=yes" | sudo tee -a /etc/samba/smb.conf
    sudo service smbd restart
    ```
12. Set a samba password for your username using `sudo smbpasswd -a $USER`
13. Allow ssh/samba connections through the firewall and then enable it:
    ```bash
    sudo ufw allow ssh
    sudo ufw allow samba
    sudo ufw enable
    ```

### Container Engine

1. Install docker and its dependencies:
    ```bash
    sudo apt install -y apt-transport-https ca-certificates curl gnupg lsb-release
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
    echo "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list
    sudo apt update
    sudo apt install -y docker-ce docker-ce-cli containerd.io
    ```
2. Add yourself to the `docker` group and then logout and login again:
    ```bash
    sudo usermod -a -G docker $USER
    exit
    ```
3. Make sure that docker works by pulling a container to run `bwa-mem2`:
    ```bash
    docker run --rm -u $(id -u):$(id -g) ghcr.io/ucladx/bwa-mem:2.2.1 bwa-mem2
    ```

### Run Scanner

[Click here](https://miso-lims.readthedocs.io/projects/runscanner) to read about Run Scanner. These are the installation steps summarized:

1. Install dependencies, open a port, and create a config file:
    ```bash
    sudo apt install -y tomcat9 maven clang cmake libjsoncpp-dev autoconf libtool build-essential
    sudo ufw allow 8080/tcp
    echo -e '<Context>\n    <Parameter name="runscanner.configFile" value="/etc/runscanner.json" override="false"/>\n</Context>' | sudo tee /var/lib/tomcat9/conf/Catalina/localhost/runscanner.xml
    ```
2. Create `/etc/runscanner.json` and use it to describe paths to run folders per sequencer. For example:
    ```json
    [
        {
            "path": "/hot/runs/miseq",
            "platformType": "ILLUMINA",
            "name": "default",
            "timeZone": "America/Los_Angeles",
            "parameters": {}
        },
        {
            "path": "/hot/runs/novaseq",
            "platformType": "ILLUMINA",
            "name": "default",
            "timeZone": "America/Los_Angeles",
            "parameters": {}
        }
    ]
    ```
3. Download, build, and deploy runscanner with support for parsing Illumina InterOp files:
    ```bash
    git clone --recursive git@github.com:ucladx/runscanner.git
    cd runscanner/runscanner-illumina
    export CXX=$(which clang++) && ./build-illumina-interop && autoreconf -i && ./configure && make && sudo make install
    cd ..
    mvn clean package
    sudo cp scanner/target/scanner-*.war /var/lib/tomcat9/webapps/runscanner.war
    ```
4. Open `localhost:8080/runscanner` in a browser to make sure the graphical interface works.
5. Now you can configure each sequencer to write runs into their respective network folder via SMB.

### User Environment

1. Install mamba (fast conda replacement that uses conda-forge as default channel):
    ```bash
    curl -L https://github.com/conda-forge/miniforge/releases/download/4.12.0-2/Mambaforge-Linux-x86_64.sh -o /tmp/mambaforge.sh
    sh /tmp/mambaforge.sh -bfp $HOME/mambaforge && rm -f /tmp/mambaforge.sh
    ```
2. Add the following lines into your `.bashrc` file, then logout and login to load it:
    ```bash
    # Add mambaforge to PATH if found, and load the base environment
    if [ -f "$HOME/mambaforge/etc/profile.d/conda.sh" ]; then
        . $HOME/mambaforge/etc/profile.d/conda.sh
        conda activate
    fi
    ```
3. Create two conda environments, one with basic bioinformatics tools, and the other for interaction with AWS/DNAnexus:
    ```bash
    mamba create -n bio -c bioconda htslib samtools bcftools bedtools ucsc-liftover
    mamba create -n aws -c bioconda dxpy aws-okta-keyman awscli
    ```
