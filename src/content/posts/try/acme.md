---
title: acme+cloudflarr获取免费证书配置教程
published: 2024-09-09
description: "acme+cloudflarr获取免费证书"
image: ""
tags: ["Linux"]
category: "指南"
draft: false
---

## 安装 acme

安装很简单, 一个命令:

```bash
curl https://get.acme.sh | sh -s email=my@example.com
# 下载失败需要配置代理
export https_proxy:http://127.0.0.1:7890
```

普通用户和 root 用户都可以安装使用。

在用户目录的`.bashrc`末尾创建一个 alias

```bash
alias acme.sh=~/.acme.sh/acme.sh
```

重新进入终端或者`source .bashrc`

## 配置环境

> Cloudflare Domain API 提供两种自动颁发证书的方法\
> （a） 创建具有特定权限的限制性 API 令牌;\
> （b） 使用与您的 Cloudflare 帐户关联的全局 API 密钥，该密钥具有所有权限。

我们选择创建用户 API 令牌方式

进入 cloudflare[用户令牌](https://dash.cloudflare.com/profile/api-tokens)页面后，点击蓝色按钮**创建令牌**

![](https://cdn.la02.cc/pichub/2024/09/09/1725828818.png)

使用 `编辑区域DNS` 模板

![](https://cdn.la02.cc/pichub/2024/09/09/1725828873.png)

选择一个区域资源点击继续之后获取令牌

![](https://cdn.la02.cc/pichub/2024/09/09/1725828948.png)

返回 linux 终端配置环境变量

```bash
# 修改为刚生成的CF令牌
export CF_Token="Y_jpG9AnfQmuX5Ss9M_qaNab6SQwme3HWXNDzRWs"
```

## 创建证书

```bash
# 创建一个泛域名证书
acme.sh --issue --dns dns_cf -d xx.com -d *.xx.com
```

这里的 cf_token 会被记录，将来你在使用 cf dns api 的时候, 就不需要再次指定了。

前面证书生成以后, 接下来需要把证书 copy 到真正需要用它的地方。默认生成的证书都放在安装目录下: `~/.acme.sh/`, 请不要直接使用此目录下的文件。

正确的使用方法是使用 --install-cert 命令,并指定目标位置, 然后证书文件会被 copy 到相应的位置, 例如:

如果你是 nginx 的话可以使用下面的命令配置证书

```bash
acme.sh --install-cert -d xxx.com \
--key-file       /path/to/nginx/key.pem  \
--fullchain-file /path/to/nginx/cert.pem \
--reloadcmd     "service nginx force-reload"
```

也可以手动复制内容修改

(一个小提醒, 这里用的是 service nginx force-reload, 不是 service nginx reload, 据测试, reload 并不会重新加载证书, 所以用的 force-reload)

在安装 acme.sh 之后会默认生成一个定时任务

```bash
18 7 * * * "/home/username/.acme.sh"/acme.sh --cron --home "/home/username/.acme.sh" > /dev/null
```

作用是每天早上 7 点 18 分执行位于`/home/username/.acme.sh`目录下的 acme.sh 脚本，以更新 SSL 证书，并且所有的输出信息都会被丢弃，不会在系统日志中留下记录。

> [源地址](https://github.com/acmesh-official/acme.sh/wiki/%E8%AF%B4%E6%98%8E)
