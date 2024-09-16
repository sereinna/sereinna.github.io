---
title: 简易指南列表
published: 2024-09-13
description: "Linux系统，WSL系统配置"
image: ""
tags: [Linux]
category: "简易指南"
draft: false
lang: ""
---

[Linux](#linux)\
[WSL](#wsl)

## Linux

[初始化](#初始化配置)\
[安装 docker](#安装-docker)

### 初始化配置

[添加用户](#添加用户)\
[设置主机名](#设置主机名)\
[更新软件包](#更新软件包)

#### 添加用户

创建普通 bash 用户 [**可选**]

```bash
useradd -m -s /bin/bash 新用户名
passwd 新用户名
```

添加新用户拥有管理员权限

修改 sudoers 文件，在 root 那一行下面添加用户

```bash
sudo sed -i '/^root/ a\用户名\tALL=(ALL:ALL) ALL' /etc/sudoers
#例如
sudo sed -i '/^root/ a\ubuntu\tALL=(ALL:ALL) ALL' /etc/sudoers
```

#### 设置主机名

```bash
sudo hostnamectl set-hostname 主机名
# 刷新终端
bash
```

#### 更新软件包

<details>

<summary>24.04 以下</summary>

```bash
sudo sh -c 'cat <<EOF > /etc/apt/sources.list
# 默认注释了源码镜像以提高 apt update 速度，如有需要可自行取消注释
deb https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs) main restricted universe multiverse
# deb-src https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs) main restricted universe multiverse
deb https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-updates main restricted universe multiverse
# deb-src https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-updates main restricted universe multiverse
deb https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-backports main restricted universe multiverse
# deb-src https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-backports main restricted universe multiverse

# 以下安全更新软件源包含了官方源与镜像站配置，如有需要可自行修改注释切换
# deb https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-security main restricted universe multiverse
# # deb-src https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-security main restricted universe multiverse

deb http://security.ubuntu.com/ubuntu/ $(lsb_release -cs)-security main restricted universe multiverse
# deb-src http://security.ubuntu.com/ubuntu/ $(lsb_release -cs)-security main restricted universe multiverse

# 预发布软件源，不建议启用
# deb https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-proposed main restricted universe multiverse
# # deb-src https://mirror.nju.edu.cn/ubuntu/ $(lsb_release -cs)-proposed main restricted universe multiverse
EOF'
```

</details>

<details>

<summary>24.04 及以上</summary>

```bash
sudo sh -c 'cat <<'EOF' > /etc/apt/sources.list.d/ubuntu.sources
Types: deb
URIs: https://mirror.nju.edu.cn/ubuntu
Suites: noble noble-updates noble-backports
Components: main restricted universe multiverse
Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg

# 默认注释了源码镜像以提高 apt update 速度，如有需要可自行取消注释
# Types: deb-src
# URIs: https://mirror.nju.edu.cn/ubuntu
# Suites: noble noble-updates noble-backports
# Components: main restricted universe multiverse
# Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg

# 以下安全更新软件源包含了官方源与镜像站配置，如有需要可自行修改注释切换
# Types: deb
# URIs: https://mirror.nju.edu.cn/ubuntu
# Suites: noble-security
# Components: main restricted universe multiverse
# Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg

# # Types: deb-src
# # URIs: https://mirror.nju.edu.cn/ubuntu
# # Suites: noble-security
# # Components: main restricted universe multiverse
# # Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg

Types: deb
URIs: http://security.ubuntu.com/ubuntu/
Suites: noble-security
Components: main restricted universe multiverse
Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg

# Types: deb-src
# URIs: http://security.ubuntu.com/ubuntu/
# Suites: noble-security
# Components: main restricted universe multiverse
# Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg

# 预发布软件源，不建议启用
# Types: deb
# URIs: https://mirror.nju.edu.cn/ubuntu
# Suites: noble-proposed
# Components: main restricted universe multiverse
# Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg

# # Types: deb-src
# # URIs: https://mirror.nju.edu.cn/ubuntu
# # Suites: noble-proposed
# # Components: main restricted universe multiverse
# # Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg
EOF'
```

</details>

### 安装 docker

添加阿里云密钥

```bash
curl -fsSL http://mirrors.aliyun.com/docker-ce/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
```

载入软件源

```bash
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] http://mirrors.aliyun.com/docker-ce/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
```

更新软件源

```bash
sudo apt update
```

安装 docker

```bash
sudo apt install docker-ce
```

配置 docker 用户

```bash
sudo usermod -aG docker $USER
# 更新用户组
newgrp docker
```

设置镜像加速

```bash
sudo vim /etc/docker/daemon.json
```

```json
{
  "registry-mirrors": [
    "https://hub.uuuadc.top",
    "https://docker.anyhub.us.kg",
    "https://dockerhub.jobcher.com",
    "https://dockerhub.icu",
    "https://docker.ckyl.me",
    "https://docker.awsl9527.cn"
  ]
}
```

## WSL

[安装 WSL](#安装-wsl-系统)\
[安装 Arch 系统](#使用-wsl-ubuntu-系统配置-arch-系统)

### 安装 WSL 系统

进入程序的启用与关闭，开启`适用于Linux的Windows子系统`和`虚拟机平台`
![图 0](https://cdn.la02.cc/pichub/2024/09/14/1726274353.png)  
![图 1](https://cdn.la02.cc/pichub/2024/09/14/1726274363.png)

启用后**重启系统**

更新 wsl 系统

```bash
wsl --update
```

设置 WSL2 系统

```bash
wsl --set-default-version 2
```

安装 ubuntu 系统

```bash
# 查看线上系统
wsl --list --online
"
使用 'wsl.exe --install <Distro>' 安装。

NAME                            FRIENDLY NAME
Ubuntu                          Ubuntu
Debian                          Debian GNU/Linux
kali-linux                      Kali Linux Rolling
Ubuntu-18.04                    Ubuntu 18.04 LTS
Ubuntu-20.04                    Ubuntu 20.04 LTS
Ubuntu-22.04                    Ubuntu 22.04 LTS
Ubuntu-24.04                    Ubuntu 24.04 LTS
OracleLinux_7_9                 Oracle Linux 7.9
OracleLinux_8_7                 Oracle Linux 8.7
OracleLinux_9_1                 Oracle Linux 9.1
openSUSE-Leap-15.6              openSUSE Leap 15.6
SUSE-Linux-Enterprise-15-SP5    SUSE Linux Enterprise 15 SP5
SUSE-Linux-Enterprise-15-SP6    SUSE Linux Enterprise 15 SP6
openSUSE-Tumbleweed             openSUSE Tumbleweed
"
wsl --install '<Distro>'
# 例如 wsl --install Ubuntu 使用前面的NAME的系统名称
```

安装 Ubuntu 需要配置用户名和密码，之后使用 `wsl.exe` 可以进入默认系统

```bash
wsl --list # 查看当前系统
wsl -d 系统名称 # 进入选定系统
```

![图 2](https://cdn.la02.cc/pichub/2024/09/14/1726274852.png)

备份系统

```bash
wsl --export [系统名称] [.tar文件路径]
# 例如
wsl --export u2404 C:/u2404.tar
```

安装备份系统

```bash
wsl --import [自定系统名] [系统位置] [系统文件目录]
# 例如
wsl --import u2404-1 C:/wsl/u2404/ C:/u2404.tar
```

:::tip[注意]
使用导入方式默认使用 root 用户\
需要修改/etc/wsl.conf 文件\
在末尾添加

```ini
[user]
default=username
```

:::

### 使用 WSL ubuntu 系统配置 Arch 系统

进入任意一个 ubuntu 系统

```bash
# 更新系统
sudo apt update
sudo apt upgrade
# 安装软件包
sudo apt install libarchive-tools
```

下载 Arch 最新系统

```bash
wget https://mirrors.bfsu.edu.cn/archlinux/iso/latest/archlinux-bootstrap-x86_64.tar.zst
sudo bsdtar xpf archlinux-bootstrap-x86_64.tar.zst
sudo bsdtar cpf archlinux.tar -C root.x86_64/ .
```

移动导入 配置镜像源

```bash
cp archlinux.tar /mnt/c/
wsl --import Arch C:/Linux/WSL/arch C:/archlinux.tar
echo 'Server = https://mirrors.nju.edu.cn/archlinux/$repo/os/$arch' > /etc/pacman.d/mirrorlist
```

设置密钥更新软件包

```bash
pacman-key --init && pacman-key --populate
pacman -Sy archlinux-keyring && pacman -Su
```
