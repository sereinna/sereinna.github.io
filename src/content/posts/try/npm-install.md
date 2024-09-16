---
title: 使用Nginx Proxy Manager配置反向代理
published: 2024-09-11
description: "反向代理安装配置教程"
image: "https://cdn.la02.cc/pichub/2024/09/11/1726017896.png"
tags: [Docker, Web项目, Linux]
category: 指南
draft: false
lang: ""
---

## 前言

:::tip
在阅读本文时，如果您有任何疑问或需要帮助，请通过以下方式联系我们：\
发送相关问题至 admin@la02.club，您的反馈将帮助博主完善文章。
:::

## 安装 docker

使用阿里云[安装 docker](/posts/server-setup/#安装-docker)，并设置加速源。

## 安装 Nginx Proxy Manager

### 首先，我们需要创建一个目录来存放 Nginx Proxy Manager 的相关配置。

```bash
mkdir ~/docker/npm & cd ~/docker/npm
```

### 创建 docker-compose.yml

使用以下命令创建并编辑`docker-compose.yml`文件：

```bash
vim docker-compose.yml
```

在文件中添加以下配置：

```yml
services:
  app:
    image: "chishin/nginx-proxy-manager-zh:release"
    restart: unless-stopped
    ports:
      # These ports are in format <host-port>:<container-port>
      - "80:80" # Public HTTP Port
      - "443:443" # Public HTTPS Port
      - "81:81" # Admin Web Port
      # Add any other Stream port you want to expose
      # - '21:21' # FTP

    # Uncomment the next line if you uncomment anything in the section
    # environment:
    # Uncomment this if you want to change the location of
    # the SQLite DB file within the container
    # DB_SQLITE_FILE: "/data/database.sqlite"

    # Uncomment this if IPv6 is not enabled on your host
    # DISABLE_IPV6: 'true'

    volumes:
      - ./data:/data
      - ./letsencrypt:/etc/letsencrypt
```

:::tip[提示]
开放 80 443 端口\
Web 站点在 81 端口\
web 站点的默认用户名密码：

- 用户名：admin@example.com
- 密码：changeme
  :::

### 启动服务

```bash
docker compose up -d
```

### 登录并配置反向代理

访问管理界面，地址为 ip:81。

:::tip[提示]
在终端输入 `curl ip.sb` 获取公网 ip

使用默认账号密码登录之后要设置新邮箱和密码
:::

#### 配置 DNS 域名解析

在域名解析服务中添加以下记录：

- `@` 记录指向公网 IP 地址
- `npm` 记录指向公网 IP 地址

![picture 0](https://cdn.la02.cc/pichub/2024/09/11/1726016619.png)

cloudflare 配置为例，添加 DNS 解析，配置`A记录类型`，`@主机记录`，`IP地址`，关闭 CF 代理。

`npm`同理，修改`@`为对应的名称

#### 创建反向代理并获取免费证书

在 NPM 仪表盘页面选择“代理服务”并点击“添加代理服务”。

![picture 2](https://cdn.la02.cc/pichub/2024/09/11/1726017205.png)

填写域名（例如 npm.xxx.com），选择默认协议 HTTP，设置目标 IP 地址和端口。

![picture 8](https://cdn.la02.cc/pichub/2024/09/11/1726019856.png)

配置 SSL 证书，选择“申请新的 SSL 证书”，勾选“开启 SSL”，同意条款后保存。

![picture 3](https://cdn.la02.cc/pichub/2024/09/11/1726017487.png)

等待 SSL 证书获取完成后，即可通过配置的域名访问服务。

![picture 4](https://cdn.la02.cc/pichub/2024/09/11/1726017640.png)

### 泛域名证书

#### 生成泛域名证书

请参考[泛域名证书生成教程](/posts/acme/)

#### 上传泛域名证书

点击`SSL证书` --> `添加SSL证书` --> `上传证书`

![picture 6](https://cdn.la02.cc/pichub/2024/09/11/1726018436.png)

有三个必填项

![picture 7](https://cdn.la02.cc/pichub/2024/09/11/1726018807.png)

|  名字\*  |   密钥\*   |   证书\*   |   中间证书    |
| :------: | :--------: | :--------: | :-----------: |
| \*.xx.xx | domain.key | domain.cer | fullchain.cer |

最后编辑代理服务的 SSL 证书来源，配置上传的证书文件即可。
