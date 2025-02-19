---
title: minecraft服务器搭建
published: 2024-09-18
description: "一个基于linux的1.21的minecraft服务器的一键搭建流程"
image: "./120709207_p0.jpg"
tags: [minecraft,linux]
category: "指南"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 mc服务器 `minecraft`
:::

## 命令行启动代码部分
```bash 

#!/bin/bash

# echo "xxx"

# 省略密码输入和验证部分
echo "密码验证已跳过"

if [ ! -e "./.flag" ]; then
    read -p "请输入您的系统版本（ubuntu（默认）/syno（群晖））：" sys_version

    if [ "$sys_version" == "syno" ]; then
        # sudo opkg update
        # sudo opkg install wget
        # sudo opkg install net-tools
        wget https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.tar.gz
        tar -zxvf jdk-21_linux-x64_bin.tar.gz
        sudo mv jdk-21.0.4 /opt/

        profile_grep=$(cat /etc/profile | grep JAVA)
        if [ "$profile_grep" != "" ]; then
            echo "环境变量已存在，跳过写入"
        else
            echo 'export JAVA_HOME="/opt/jdk-21.0.4"' | sudo tee -a /etc/profile
            echo 'export PATH="$JAVA_HOME/bin:$PATH"' | sudo tee -a /etc/profile
        fi

        sudo touch ./.flag
        read -p "基础环境搭建成功，回车将重启，请重启之后再进入路径并执行此脚本。"
        sudo init 6
    else
        sudo apt update
        sudo apt install openjdk-21-jdk
        sudo apt install wget
        sudo apt install net-tools
        sudo ufw allow 25565
    fi
fi

read -p "请选择您要安装的服务端(官服(o)/fabric(f)/mohist(m，支持forge的第三方服务端))：" type
read -p "请输入您要安装的minecraft版本：" version

if [ "$type" == "o" ]; then
    # 从mcversions.net获取对应mc版本的官方服务端下载地址
    url=$(curl https://mcversions.net/download/$version | sed 's/<\//
</g' | grep "Download Server Jar" | grep -oP 'href="\K[^"]+')
    wget $url
elif [ "$type" == "m" ]; then
    wget -O "server.jar" https://mohistmc.com/api/$version/latest/download
elif [ "$type" == "f" ]; then
    echo "fabric默认版本0.15.11"
    wget -O "server.jar" https://meta.fabricmc.net/v2/versions/loader/$version/0.15.11/1.0.0/server/jar
else
    echo "输入有误！"
    exit
fi

echo "注:部分版本可能需要您主动将目录下的eula文件改成true，请手动修改后运行run.sh运行服务器"
read -p "请输入您想要设置的内存大小（例：1024M）：" memory
touch run.sh
echo "java -Xmx$memory -Xms$memory -jar server.jar nogui" > run.sh
chmod u+x run.sh
echo "搭建成功，请手动开放防火墙并运行./run.sh，脚本将自动清理"
rm -rf ./.flag
rm -rf jdk-21_linux-x64_bin.tar.gz
rm -rf ./minecraft-init.sh

```