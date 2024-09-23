---
title: 阿里云服务器版本设置
published: 2024-09-23
description: "如何在服务器上挂在自己的阿里云盘"
image: "./117157990_p0_master1200.jpg"
tags: [网盘,linux,存储]
category: "指南"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 阿里云账号 `会员账号`
:::

## 源代码网页
https://github.com/tickstep/aliyunpan



## 正常安装步骤

### github安装
```bash
wget https://github.com/tickstep/aliyunpan/releases/download/v0.3.2/aliyunpan-v0.3.2-linux-amd64.zip

unzip aliyunpan-v0.3.2-linux-amd64.zip

cd aliyunpan-v0.3.2-linux-amd64

./aliyunpan
```
### apt安装
```bash
sudo curl -fsSL http://file.tickstep.com/apt/pgp | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/tickstep-packages-archive-keyring.gpg > /dev/null && echo "deb [signed-by=/etc/apt/trusted.gpg.d/tickstep-packages-archive-keyring.gpg arch=amd64,arm64] http://file.tickstep.com/apt aliyunpan main" | sudo tee /etc/apt/sources.list.d/tickstep-aliyunpan.list > /dev/null && sudo apt-get update && sudo apt-get install -y aliyunpan
```

## 使用相关

```bash
cd aliyunpan-v0.3.2-linux-amd64
./aliyunpan
aliyunpan > help
```


| 命令           | 描述                           |  
| -------------- | --------------------------------- |  
| `album, abm`  | 相簿(Beta)                     |  
| `cd`          | 切换工作目录                   |  
| `download, d` | 下载文件/目录                  |  
| `ls, l, ll`   | 列出目录                       |  
| `mkdir`       | 创建目录                       |  
| `mv`          | 移动文件/目录                  |  
| `pwd`         | 输出工作目录                   |  
| `recycle`     | 回收站                         |  
| `rename`      | 重命名文件                     |  
| `rm`          | 删除文件/目录                  |  
| `share`       | 分享文件/目录                  |  
| `sync`        | 同步备份功能                   |  
| `upload, u`   | 上传文件/目录                  |