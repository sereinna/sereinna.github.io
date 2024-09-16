---
title: 使用nginx配置一个HTTPSWeb站点
published: 2024-09-12
description: "使用nginx配置一个HTTPSWeb站点"
image: ""
tags: [Linux, Docker, Web项目]
category: "指南"
draft: false
lang: ""
---

:::tip[注意]
本文档中的代码顺序可能与视频教程中的顺序不完全一致，请根据视频中的操作步骤仔细核对代码的作用。
:::

## 创建目录并进入

首先，我们需要创建一个目录结构来存放 Nginx 相关的配置文件和证书。

`mkdir -p docker/nginx/{ssl,html,conf.d} & cd docker/nginx/`

## 创建 docker-compose.yml 文件

接下来，我们将创建一个 docker-compose.yml 文件来定义 Nginx 服务。

`vim compose.yml`

在 `compose.yml`文件中写入以下配置：

```yml
services:
  nginx:
    image: nginx:latest
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./ssl:/etc/nginx/ssl
      - ./html:/usr/share/nginx/html
      - ./conf.d/:/etc/nginx/conf.d/
    restart: always
```

### 创建 `Nginx` 配置文件和 `HTML` 首页文件

```bash
touch conf.d/1.conf
echo "Welcome to Nginx!" > html/index.html
```

### 为 nginx 安装证书

```bash
acme.sh --install-cert -d example.com \
--key-file       /home/ubuntu/docker/nginx/ssl/server.key  \
--fullchain-file /home/ubuntu/docker/nginx/ssl/server.cer
```

## 编写 nginx 配置文件

```conf
server {
    listen       80;
    listen       443 ssl;
    listen  [::]:80;
    listen  [::]:443 ssl;
    server_name  example.com;

    ssl_certificate /etc/nginx/ssl/server.cer;
    ssl_certificate_key /etc/nginx/ssl/server.key;

    ssl_session_timeout 1d;
    ssl_session_cache shared:MozSSL:10m;
    ssl_session_tickets off;

    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers 'ECDHE-ECDSA-AES128-GCM-SHA256:ECDHE-RSA-AES128-GCM-SHA256:ECDHE-ECDSA-AES256-GCM-SHA384:ECDHE-RSA-AES256-GCM-SHA384:DHE-RSA-AES128-GCM-SHA256:DHE-RSA-AES256-GCM-SHA384';
    ssl_prefer_server_ciphers off;

    add_header Strict-Transport-Security "max-age=31536000" always;

    location / {
        root   /usr/share/nginx/html;
        index  index.html index.htm;
    }

    error_page   500 502 503 504  /50x.html;
    location = /50x.html {
        root   /usr/share/nginx/html;
    }

}
```

## 启动测试

最后，使用以下命令启动 Nginx 服务进行测试。如果您希望 Nginx 在后台运行，可以使用`-d` 选项。

```bash
docker compose up
```

或者后台运行：

```bash
docker compose up -d
```

请确保在执行以上命令前，您已经正确配置了所有必要的文件和目录。
