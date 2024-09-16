---
title: Ubuntu安装Teamspeak服务,自定义端口
published: 2024-09-13
description: "Ubuntu安装Teamspeak服务"
image: ""
tags: [Linux]
category: "指南"
draft: false
lang: ""
---

更新软件包，安装相关软件

```bash
sudo apt update
sudo apt upgrade
sudo apt install bzip2 inetutils-telnet
```

新建一个系统用户管理服务

```bash
sudo useradd -r -m -U -d /opt/teamspeak -s /bin/bash teamspeak
sudo su - teamspeak
```

下载软件包

```bash
wget https://files.teamspeak-services.com/releases/server/3.13.7/teamspeak3-server_linux_amd64-3.13.7.tar.bz2
```

解压重命名，进入目录

```bash
tar xjf teamspeak3-server_linux_amd64-3.13.7.tar.bz2 && \
mv teamspeak3-server_linux_amd64 server && cd server
```

同意许可

```bash
touch .ts3server_license_accepted
```

:::warning
第一次启动要保存好日志中的管理员信息\
保存完毕后可以使用`Ctrl+C`退出，会在后台运行
:::

```bash
./ts3server_startscript.sh start
```

进入控制台，修改配置

:::tip
第一次启动的日志保存着管理员信息\
`loginname= "serveradmin", password= "UuHIMeKW"`
:::

```bash
telnet 127.0.0.1 10011
login serveradmin password
```

查看服务器，免费版本应该只有一个服务器

```bash
serverlist
```

![图 0](https://cdn.la02.cc/pichub/2024/09/14/1726280971.png)

修改 server_id=1 的服务器信息

```bash
use sid=1 -virtual
serveredit virtualserver_port=19987
use sid=0
```

修改文件传输端口

```bash
instanceedit serverinstance_filetransfer_port=50033
```

重启

```bash
serverstop sid=1
serverstart sid=1
```

查看修改

```bash
serverlist
instanceinfo
```

详细过程

![图 1](https://cdn.la02.cc/pichub/2024/09/14/1726281264.png)

修改完成后，键入 `quit` 退出

退出 telnet 之后再停止服务并返回 root 用户

```bash
./ts3server_startscript.sh stop
exit
```

创建守护进程文件

```bash
vim /etc/systemd/system/teamspeak.service
```

编辑文件

```ini
[Unit]
Description=TeamSpeak3 Server
Wants=network-online.target
After=syslog.target network.target

[Service]
WorkingDirectory=/opt/teamspeak/server/
User=teamspeak
Group=teamspeak
Type=forking
ExecStart=/opt/teamspeak/server/ts3server_startscript.sh start
ExecStop=/opt/teamspeak/server/ts3server_startscript.sh stop
ExecReload=/opt/teamspeak/server/ts3server_startscript.sh restart
PIDFile=/opt/teamspeak/server/ts3server.pid

[Install]
WantedBy=multi-user.target
```

服务管理

```bash
systemctl daemon-reload
systemctl start teamspeak
systemctl status teamspeak
systemctl enable teamspeak
systemctl stop teamspeak
```
