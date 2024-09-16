---
title: ubuntu创建Minecraft服务器
published: 2024-08-27
description: "这个指南解释如何在 Ubuntu 20.04 上如何搭建我的世界服务器。我们将会使用 Systemd 来运行我的世界服务器以及mcrcon工具来连接运行的实例。"
image: ""
tags: ["服务器", "Linux"]
category: "指南"
draft: false
---

Minecraft 一直是最流行的游戏之一。它是一个沙盒视频游戏，用户可以体验无限的世界，并且可以构建不同的结构，从简单的房子到高耸的摩天大楼。这个指南解释如何在 Ubuntu 上如何搭建我的世界服务器。我们将会使用 Systemd 来运行我的世界服务器以及 mcrcon 工具来连接运行的实例。我们也将向你展示如何创建一个计划任务，执行常规的服务器备份。

## 一、前提

根据 Minecraft 官方网站，4GB RAM 内存是最基本的配置。

安装必要的软件包来构建 mcrcon 工具：

```bash
sudo apt update
sudo apt install git build-essential
```

## 二、安装运行环境 Java mcrcon

### JAVA

Minecraft 需要 Java 8 或者更高版本。 我们将会安装 Java 的 JRE 版本。这个版本更适合服务器应用，因为它有更少的依赖，并且使用更少的系统资源。

我们使用压缩包安装

进入 java 构建网站[**adoptium**](https://adoptium.net/zh-CN/temurin/releases/?package=jre&os=linux&arch=x64&version=21)，当前链接已设定 jre21、、linux、x64，可手动修改

```bash
# 下载完成后使用scp上传到服务器
scp ./OpenJDK*.tar.gz username@host:~/
# 移动到/usr/local/
sudo mv OpenJDK*.tar.gz /usr/local/ && cd /usr/local
# 解压文件
sudo tar xzf OpenJDK*.tar.gz
# ls -l 查看解压目录
# 创建软连接到/usr/local/bin目录，并重命名放置冲突
sudo ln -s /usr/local/jdk-21.0.4+7-jre/bin/java /usr/local/bin/jre21
# 查看软连接
ls -l /usr/local/bin/jre21
```

### mcrcon

RCON 是一个协议，它允许你连接到 Minecraft 服务器，并且执行命令。[mcron](https://github.com/tiiffi/mcrcon)是一个 RCON 客户端，使用 C 语言编写而成。

我们将会从 Github 下载最新的源代码，并且构建 `mcrcon`二进制文件。

从 Github 克隆 `Tiiffi/mcrcon`源到 `~/tools/mcron`目录：

```bash
# 无法下载可以尝试使用代理 git config --global https.proxy http://127.0.0.1:7890
# 自行切换代理地址，或使用其他github加速方式
sudo git clone https://github.com/Tiiffi/mcrcon.git /usr/local/mcrcon
cd /usr/local/mcrcon
make
```

一旦完成，验证 `mcrcon`编译成功，打印它的版本：

```bash
./mcrcon -v
```

```输出
mcrcon 0.7.2 (built: Feb 22 2024 01:39:23) - https://github.com/Tiiffi/mcrcon
Bug reports:
        tiiffi+mcrcon at gmail
        https://github.com/Tiiffi/mcrcon/issues/
```

## 三、创建 Minecraft 用户

因为安全原因， Minecraft 不应该在 root 用户下运行。我们将会创建一个新的系统用户和用户组，用户主目录 `/opt/minecraft`。这个用户有最小权限，来运行 Minecraft 服务器：

```bash
sudo useradd -r -m -U -d /opt/minecraft -s /bin/bash minecraft
```

我们不会为这个用户设置密码。这样，这个用户将不能通过 SSH 登录。想要修改 `minecraft`用户，你将需要使用 root 登录服务器，或者其他有 sudo 权限的用户。

## 四、在 Ubuntu 上安装 Minecraft

在开始安装过程之前，切换到 minecraft 用户：

```bash
sudo su - minecraft
```

复制
运行下面的命令在用户主目录下创建三个新的目录：

```bash
mkdir -p ~/{backups,tools,server}
```

- `backups`目录将会存储 服务器备份。你可以同步这个目录到你的远程备份服务器。
- `tools` 目录将会托管 `mcrcon`客户端和备份脚本。
- `server`目录将会包含实际的 Minecraft 服务器和它的数据。

### 4.1 下载 Minecraft 服务器

有一些 Minecraft 服务器 mods 例如 [**Craftbukkit**](https://getbukkit.org/download/craftbukkit) 或者 [**Spigot**](https://www.spigotmc.org/) ，允许你在你的服务器上添加特性（插件）以及定制，以及调整服务器设置。

在这个指南中，我们将会安装最新的 Mojang 官方 vanilla 我的世界服务器。同样的指令，同样适合于其他的服务器 mods。

去 [Minecraft 下载页面](https://minecraft.net/en-us/download/server/) 下载最新的 Minecraft 服务器 Java 压缩包（JAR）。在写作的时候，最新的版本是 1.21.1。我们使用[**Paper**](https://papermc.io/downloads/paper)核心创建服务器。

使用 `wget`下载 jar 文件到~/server 目录。

```bash
wget https://api.papermc.io/v2/projects/paper/versions/1.21.1/builds/52/downloads/paper-1.21.1-52.jar -OP ~/server/server.jar
```

### 4.2 配置 Minecraft 服务器

一旦下载完成，切换到 `~/server`目录，并且启动 Minecraft 服务器：

```shell
cd ~/server
jre21 -Xmx2048M -Xms1024M -jar server.jar nogui
```

第一次启动的时候，服务器执行一些操作，创建 `server.properties`和 `eula.txt`文件，并且停止。

```输出
[17:35:14] [main/ERROR]: Failed to load properties from file: server.properties
[17:35:15] [main/WARN]: Failed to load eula.txt
[17:35:15] [main/INFO]: You need to agree to the EULA in order to run the server. Go to eula.txt for more info.
```

想要运行服务器，你需要同意 `Minecraft EULA`，就像上面输出所有提示的。

打开 `eula.txt`文件，并且修改 `eula=false` 为 `eula=true`:

```bash
vim ~/server/eula.txt
```

```txt
# false修改为true
eula=true
```

关闭并且保存文件。

下一步，打开 `server.properties`文件，并且启动 `rcon`协议，并且设置 rcon 密码：

```bash
vim ~/server/server.properties
```

定位到下面的行，并且更新它们的值，就像下面显示的一样：

```properties
rcon.port=25575
rcon.password=strong-password
enable-rcon=true
```

不要忘记将 `strong-password`修改为一些更安全的密码。如果你不想从远程位置访问 Minecraft 服务器，确保 rcon 端口被防火墙所阻塞。
复制
在这里，你可以调整服务器的默认属性。想要获取更多关于服务器设置的信息，浏览[minecraft.service](https://minecraft.gamepedia.com/Server.properties)页面。

## 五、创建 Systemd 单元文件

与手动启动 Minecraft 服务器相比，我们将会创建一个 Systemd 单元文件，并且将 Minecraft 当作服务来运行。

通过 `exit`来切换回你的 sudo 用户。

打开你的文本编辑器，并且在 `/etc/systemd/system/`目录下创建一个名为 `mcpaper.service`的文件：

```bash
sudo vim /etc/systemd/system/mcpaper.service
```

粘贴下面的配置：

```ini
[Unit]
Description=Minecraft Server
After=network.target

[Service]
User=minecraft
Nice=1
KillMode=none
SuccessExitStatus=0 1
ProtectHome=true
ProtectSystem=full
PrivateDevices=true
NoNewPrivileges=true
WorkingDirectory=/opt/minecraft/server
ExecStart=jre21 -Xmx2048M -Xms1024M -jar server.jar nogui
ExecStop=mcrcon -H 127.0.0.1 -P 25575 -p strong-password stop

[Install]
WantedBy=multi-user.target
```

根据你的服务器资源来调整 `Xmx`和 `Xms`标志。`Xmx`标志定义 Java 虚拟机的最大申请内存。而 `Xms`定义了初始申请内存。当前，确保你使用了正确的 `rcon`端口和密码。

保存文件，并且重新加载 systemd 管理配置：

```bash
sudo systemctl daemon-reload
sudo systemctl start mcpaper
sudo systemctl status mcpaper
```

输出：

```bash
● mcpaper.service - Minecraft Server
     Loaded: loaded (/etc/systemd/system/mcpaper.service; disabled; preset: enabled)
     Active: active (running) since Tue 2024-08-27 14:23:02 CST; 1min 6s ago
   Main PID: 55281 (jre21)
      Tasks: 40 (limit: 4697)
     Memory: 1.3G ()
        CPU: 38.196s
     CGroup: /system.slice/mcpaper.service
             └─55281 jre21 -Xmx2048M -Xms1024M -jar server.jar nogui
```

最后，启动 Minecraft 服务开机启动。

```bash
sudo systemctl enable mcpaper
```

## 六、调整防火墙

Ubuntu 附带防火墙工具 UFW。如果在你的系统上启用了防火墙，你想从你的本地网络访问 Minecraft 服务器，你需要打开端口 `25565`:

```bash
sudo ufw allow 25565/tcp
```

云服务器需要去对应厂商的防火墙控制台配置规则

## 七、配置备份

在这一节，我们创建一个备份 shell 脚本和计划任务，以便自动备份 Minecraft 服务器。

切换到 `minecraft`:

```bash
sudo su - minecraft
```

打开你的文本编辑器，并且创建下面的文件：

```bash
vim /opt/minecraft/tools/backup.sh
```

粘贴下面的配置：

```sh
#!/bin/bash

function rcon {
                                # 自处需修改登录密码
  mcrcon -H 127.0.0.1 -P 25575 -p strong-password "$1"
}

rcon "save-off"
rcon "save-all"
tar -cvpzf /opt/minecraft/backups/server-$(date +%F-%H-%M).tar.gz /opt/minecraft/server
rcon "save-on"

## Delete older backups
find /opt/minecraft/backups/ -type f -mtime +7 -name '*.gz' -delete
```

保存文件，并且将脚本设置为可执行：

```bash
chmod +x /opt/minecraft/tools/backup.sh
```

下一步，创建一个定时任务，每天自动在一个固定时间运行一次。

打开你的 crontab 文件，输入：

```bash
crontab -e
```

每天 23：00 运行备份脚本，粘贴下面的文本：

```bash
0 23 * * * /opt/minecraft/tools/backup.sh
```

## 八、访问 Minecraft 终端

想要访问 Minecraft 终端，使用 `mcrcon`工具。你需要指定主机，rcon 端口，rcon 密码并且使用 `-t`（启动 `mcrcon`终端模式）：

```bash
mcrcon -H 127.0.0.1 -P 25575 -p strong-password -t
```

```输出
Logged in. Type "Q" to quit!
>
```

从远程位置访问 Minecraft 终端，确保 rcon 端口没有被阻塞。

如果你正常连接到 Minecraft 终端，不想输入一大串命令，你可以创建一个 [bash 关联](https://linuxize.com/post/how-to-create-bash-aliases/)。

## 九、总结

我们已经向你展示如何在 Ubuntu 20.04 上搭建一个 Minecraft（我的世界）服务器，并且设置每天备份。

现在你可以启动你的 Minecraft 客户端，连接到服务器，并且开始 Minecraft 冒险。

原文 ：[https://linuxize.com/post/how-to-make-minecraft-server-on-ubuntu-20-04/](https://linuxize.com/post/how-to-make-minecraft-server-on-ubuntu-20-04/)

## 十、基岩版

```bash
cd ~
mkdir Server
cd Server
sudo apt-get install wget -y
wget https://minecraft.azureedge.net/bin-linux/bedrock-server-1.20.71.01.zip
```

```bash
sudo apt-get install unzip -y
unzip bedrock-server-1.20.71.01.zip
```

随后你可以使用 ls 命令查看目录下的文件，有可能这些文件被再次纳入了一个文件夹，你可以使用 cd 命令进入此文件夹，然后执行以下任一命令启动它：

```bash
LD_LIBRARY_PATH=. ./bedrock_server
或
./bedrock_server
```

使 BDS 在后台运行编辑编辑源代码
当你退出 SSH 终端亦或者结束掉 bash 终端时，你运行的程序也会随之终止，为此，你需要安装 screen：

`sudo apt-get install screen -y`

然后使用 screen 进入 screen 终端（你可能需要再按一次 ↵ Enter）

随后再次执行刚才的命令：

```bash
LD_LIBRARY_PATH=. ./bedrock_server
或
./bedrock_server
```

这样你的 BDS 就可以在后台运行而不需要你一直保持着终端开启。
