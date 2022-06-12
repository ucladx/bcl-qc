### Purpose

A proof-of-concept server to perform primary and secondary NGS analyses with reasonable TAT and a USD $2500 budget.

### Hardware and OS

Acquired a [Dell XPS 8950](https://www.dell.com/en-us/work/shop/desktops-all-in-one-pcs/xps-desktop/spd/xps-8950-desktop) tower desktop in mid 2022 with the following specs. To support ECC memory, we could have gotten a [Dell Precision 3660](https://www.dell.com/en-us/work/shop/desktops-all-in-one-pcs/precision-3660-tower-workstation/spd/precision-3660-workstation), but that goes over budget.

- Intel Core i9-12900K (best single-thread performance; added liquid cooling option)
- 64GB DDR5-4400 Memory
- 1x 960GB Intel Optane 905P PCIe NVMe SSD (nvme2n1)
- 2x 2TB Samsung 980 Pro PCIe NVMe SSD (nvme0n1, nvme1n1)

We'll choose to install Clear Linux for the latest kernel optimized for performance on Intel CPUs.

1. On a Windows PC, download the [Clear Linux Desktop ISO](https://clearlinux.org/downloads) image.
2. Boot into BIOS on the server (press F2 during Dell logo), and ensure that:
    - `SATA Operation` is set to AHCI (disables hardware RAID)
    - `Secure Boot` is disabled (does not work with Clear Linux)
    - `AC Recovery` is set to `Last Power State`
    - `POST Behaviour` is set to `Continue on Warnings` and `Keyboard Error Detection` is disabled
    - `Virtualization`, `VT for Direct I/O`, and `Trusted Execution` are enabled
3. Follow [these instructions](https://docs.01.org/clearlinux/latest/get-started/bootable-usb.html) to create a bootable USB stick.
4. Plug in the USB stick and boot into it (press F12 during Dell logo).
5. After the live desktop image boots, click "Clear Linux OS installer" on the dock and follow instructions.
6. Use the smaller SSD (nvme2n1) as the installation media and add yourself as a sudo user.
7. Click "Advanced Options" to set the hostname, make sure other options are set to your preference, and click "Install".
8. After booting into the installed OS, start Terminal and install the following:
    ```bash
    sudo swupd update
    sudo swupd bundle-add -y openssh-server storage-utils iptables clr-network-troubleshooter
    ```
9. Create a unix group for your team (E.g. `dx` below) and make it your primary group:
    ```bash
    sudo groupadd dx
    sudo usermod -g dx $USER
    sudo groupdel $USER
    ```
10. Create a level 0 RAID using the two larger SSDs and make it writable for your team:
    ```bash
    sudo mdadm --create --verbose --level=0 --raid-devices=2 /dev/md0 /dev/nvme0n1 /dev/nvme1n1
    sudo mkfs.ext4 -F /dev/md0
    sudo blkid -o value /dev/md0 | head -n1 | xargs -I{} echo "UUID={} /hot ext4 defaults,noatime 0 2" | sudo tee -a /etc/fstab
    sudo mount /hot
    sudo chown $USER:dx /hot
    sudo chmod 775 /hot
    ```
11. Use Samba and NetBIOS to share the folder `/hot` over the network:
    ```bash
    sudo mkdir -p /etc/samba
    echo -e "[hot]\n  path=/hot\n  writeable=yes" | sudo tee -a /etc/samba/smb.conf
    sudo systemctl enable --now smb
    sudo systemctl enable --now nmb
    ```
12. Set a samba password for your username using `sudo smbpasswd -a $USER`
13. Enable the firewall after allowing OpenSSH, NetBIOS, Samba, ICMP, localhost, and all server initiated connections:
    ```bash
    sudo iptables -A INPUT -p tcp -m multiport --dports 22,139,445 -j ACCEPT
    sudo iptables -A INPUT -p icmp -j ACCEPT
    sudo iptables -A INPUT -m conntrack --ctstate ESTABLISHED,RELATED -j ACCEPT
    sudo iptables -A INPUT -i lo -j ACCEPT
    sudo iptables -A OUTPUT -o lo -j ACCEPT
    sudo iptables -P OUTPUT ACCEPT
    sudo iptables -P INPUT DROP
    sudo iptables -P FORWARD DROP
    sudo systemctl start iptables-save
    sudo systemctl enable iptables-restore
    ```
