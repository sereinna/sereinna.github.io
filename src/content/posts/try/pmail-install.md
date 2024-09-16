---
title: 使用PMail创建SSL邮局(465 995)
published: 2024-09-11
description: ""
image: ""
tags: [Linux, Web项目, Docker]
category: 指南
draft: false
lang: ""
---

## 前言

:::tip
准备工作 Cloudflare `域名`，NPM `反向代理`，SSL `证书`\
相关文章:[反向代理](/posts/npm-install/)|[SSL 证书](/posts/acme/)

需放行 `465`、`995`、`8082`端口，其中`8082`为 web 端口，可修改
:::

## 创建域名解析

添加`A记录`，`mail`为主机的域名解析，关闭 CF 代理

![图 0](https://cdn.la02.cc/pichub/2024/09/11/1726041937.png)

## 配置反向代理

进入 NPM 后台添加代理 配置基础信息和 SSL 证书，证书获取请查看[前言](#前言)的相关文章

![图 2](https://cdn.la02.cc/pichub/2024/09/11/1726042222.png)

## 创建 docker-compose.yml 文件

创建一个目录保存配置文件

`mkdir -p ~/docker/pmail & cd ~/docker/pmail`

创建 yml 文件`vim docker-compose.yml`，并创建非 root 权限的配置文件夹`mkdir config`

设置环境 docker 用户组环境变量 `echo 'export DOCKER_USER="$(id -u):$(id -g)"' >> ~/.bashrc`

刷新 `source ~/.bashrc`

编辑文件，使用南京大学镜像源下载更快，`8082`端口为 Web 页面反向代理端口，可自行修改。

```yml
services:
  pmail:
    container_name: pmail
    image: "ghcr.nju.edu.cn/jinnrry/pmail:latest"
    volumes:
      - "./config:/work/config"
    ports:
      - "465:465"
      - "8082:80"
      - "995:995"
    user: "${DOCKER_USER}"
```

启动项目`docker compose up -d`

直接访问 mail.xxx.xxx 进入配置即可

![图 3](https://cdn.la02.cc/pichub/2024/09/11/1726042966.png)

之后按要求配置`数据库`，`管理员账户`、`smtp域名`和`web域名`

紧接着进入 DNS 后台解析相关配置

最后添加 SSL 证书时选择手动配置

导入`acme`生成的证书

```bash
acme.sh --install-cert -d la02.cc \
--key-file       /home/ubuntu/docker/pmail/config/ssl/private.key  \
--fullchain-file /home/ubuntu/docker/pmail/config/ssl/public.crt
```

## 测试端口放行情况

随便进入一个端口测试网站，输入 ip 地址或域名，测试 465 和 995 端口放行情况，必须要两个开放才算成功
