---
title: 服务器初始化
published: 2024-09-06
description: "云服务器Ubuntu初始化。创建普通用户，设置主机名，配置镜像源，安装docker，SSH终端汉化"
image: "https://cdn.la02.cc/pichub/2024/09/10/1725961889.png"
tags: ["Linux", "服务器"]
category: "指南"
draft: false
---

# 云服务器服务器初始化

[设置主机名](#设置主机名)\
[更新软件包](#更新软件包)\
[设置 SSH 连接后终端汉化](#设置-ssh-连接后终端汉化)\
[安装 docker](#安装-docker)

## 创建普通 bash 用户

```bash
useradd -m -s /bin/bash 新用户名
passwd 新用户名
```

- -m：创建用户的同时创建用户的家目录。
- -s /bin/bash：指定用户登录时使用的 shell，这里指定为 bash。

### 设置 sudo 权限

```bash
visudo

# 在 sudoers 文件中找到如下类似的一行：
# root    ALL=(ALL:ALL) ALL
# 在这一行下面，添加以下内容来授予用户 sudo 权限：
新用户名    ALL=(ALL:ALL) ALL
```

## 设置主机名

```bash
sudo hostnamectl set-hostname 主机名
# 刷新终端
bash
```

## 更新软件包

### 切换软件镜像源

[镜像仓库地址](https://mirror.nju.edu.cn/mirrorz-help/ubuntu/)

![图 0](https://cdn.la02.cc/pichub/2024/09/08/1725776657.png)

选择一个镜像网站，例如我们使用清华源

![图 1](https://cdn.la02.cc/pichub/2024/09/08/1725776746.png)

> 需要注意，从 Ubuntu 24.04 开始，Ubuntu 的软件源配置文件变更为 `DEB822 格式`，路径为 `/etc/apt/sources.list.d/ubuntu.sources`。

24.04 查看下方`DEB822 格式`的配置教程。

按照教程修改`/etc/apt/sources.list.d/ubuntu.sources`配置文件即可

**快速配置需要 root 用户**

![图 2](https://cdn.la02.cc/pichub/2024/09/08/1725776946.png)

24.04 之前使用`传统格式`需要切换对应 Ubuntu 版本。复制对应版本配置并修改`/etc/apt/sources.list`

**快速配置需要 root 用户**

![图 3](https://cdn.la02.cc/pichub/2024/09/08/1725777063.png)

```bash
# 更新软件包
sudo apt update
sudo apt upgrade
```

## 设置 SSH 连接后终端汉化

### 安装并启用中文

```bash
sudo apt install language-pack-zh-hans
# 查看是否已开启中文
cat /etc/locale.gen | grep zh_CN.UTF-8\ UTF-8
# 输出结果行首有 # 则未开启
# 修改文件，启用中文
sudo sed -i '/^# zh_CN.UTF-8 UTF-8/s/^# //' /etc/locale.gen
# 查看是否修改成功，行首没有 # 为开启成功
cat /etc/locale.gen | grep zh_CN.UTF-8\ UTF-8
# 重新生成 locale
sudo locale-gen
```

### 配置环境变量，检查是否为 SSH 客户端

```bash
# 一键配置
sudo sh -c "cat >> /etc/profile << 'EOF'
# 检查是否为SSH会话
if [ -n \"\$SSH_CLIENT\" ] || [ -n \"\$SSH_TTY\" ]; then
    # 设置中文环境
    export LC_ALL=zh_CN.UTF-8
    export LANG=zh_CN.UTF-8
    export LANGUAGE=zh_CN:zh
fi
EOF"
# 或者手动复制
# 检查是否为SSH会话
if [ -n "$SSH_CLIENT" ] || [ -n "$SSH_TTY" ]; then
    # 设置中文环境
    export LC_ALL=zh_CN.UTF-8
    export LANG=zh_CN.UTF-8
    export LANGUAGE=zh_CN:zh
fi
# 查看是否在末尾正确添加
cat /etc/profile
```

重新进入终端后输入`env|grep zh_CN`查看是否正确配置

## 安装 docker

### 添加 docker 密钥

```bash
curl -fsSL http://mirrors.aliyun.com/docker-ce/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
```

### 添加阿里云 docker 软件源

```bash
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] http://mirrors.aliyun.com/docker-ce/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
```

### 更新软件源

```bash
sudo apt update
```

### 安装 docker

```bash
sudo apt install docker-ce
```

### 配置 docker 用户

```bash
sudo usermod -aG docker $USER
# 更新用户组
newgrp docker
```

### 配置 docker

### 设置镜像加速

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
