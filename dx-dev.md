### Purpose

A server with sufficient compute and storage for local development, so that we don't need to touch any resources in production.

### Hardware and OS

Acquired a [Dell Precision 5820](https://www.dell.com/en-us/work/shop/desktops-all-in-one-pcs/precision-5820-tower-workstation/spd/precision-5820-workstation) tower workstation in mid 2018 with the following specs.

- Intel Xeon W-2145 (support for ECC memory and AVX-512; decent single-thread performance)
- 200GB DDR4-2666 ECC Memory (ECC reduces odds of data corruption)
- 1x 1TB Samsung 960 EVO PCIe NVMe SSD (nvme0n1)
- 1x 2TB Micron 3400 PCIe NVMe SSD (nvme1n1)
- 2x 16TB 7200rpm SATA 3.5" HDD (sda, sdb)

We'll choose to install Ubuntu 22.04 for the newer kernel and for security updates supported till 2027.

1. On a Windows PC, download the [Ubuntu Desktop 22.04 ISO](https://mirrors.ocf.berkeley.edu/ubuntu-releases/22.04) image.
2. Boot into BIOS on the server (press F2 during Dell logo), and ensure that:
    - `SATA Operation` is set to AHCI (disables hardware RAID)
    - `Secure Boot` is enabled (works with Ubuntu for now)
    - `AC Recovery` is set to `Last Power State`
    - `POST Behaviour` is set to `Continue on Warnings` and `Keyboard Error Detection` is disabled
    - `Virtualization`, `VT for Direct I/O`, and `Trusted Execution` are enabled
3. Follow [these instructions](https://ubuntu.com/tutorials/create-a-usb-stick-on-windows) to create a bootable USB stick.
4. Plug in the USB stick and boot into it (press F12 during Dell logo).
5. Choose to `Install Ubuntu` and select `Minimal installation` when prompted.
6. Under `Installation Type` choose `Something else`, and partition the first SSD (nvme0n1) as follows:
    - Two partitions, 128MB for EFI and remainder `ext4` mounted as `/`
    - Make sure the boot loader is installed on this SSD (nvme0n1)
7. Follow on-screen instructions for remaining steps (includes setting a hostname) until reboot.
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

### Container Engine

1. Install docker and its dependencies:
    ```bash
    sudo apt install -y apt-transport-https
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
3. Make sure that docker works by pulling and running `bcl-convert` to display its usage manual:
    ```bash
    docker run --rm -u $(id -u):$(id -g) ghcr.io/ucladx/bcl-convert:3.10.5 bcl-convert --help
    ```

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
    mamba create -yn bio -c bioconda htslib samtools bcftools bedtools ucsc-liftover
    mamba create -yn aws -c bioconda dxpy aws-okta-keyman awscli
    ```
