---
title: 设置网页为google可搜索
published: 2024-09-30
description: "如何把自己的网页设置为google站点地图可见"
image: "./122670348_p0.png"
tags: [网页]
category: "指南"
draft: false
lang: ""
---

## 前言  

:::tip  
准备工作 xml文件 `sitemap.xml`，HTML验证文件 `googleexxxxxx.html`，robots`robots.txt`  

:::  

## HTML 权限验证文件  

**步骤 1:** 下载验证文件  
将您的网站添加到 [Google Search Console](https://search.google.com/search-console) 后，您应该会在“推荐的验证方法”选项卡中看到一个用于下载 HTML 验证文件的选项：单击“下载文件”旁边的按钮，将此文件保存。

**步骤 2:** 上传文件  
然后把文件上传到 `public/googleexxxxxx.html`。文件位于网站的该文件夹中后，返回到 Google Search Console。  

**步骤 3:** 单击 Google Search Console 中的验证按钮  
将文件上传到您的网站后，请返回 Google Search Console 并单击“验证”按钮以完成该过程。Google Search Console 将在服务器上找到该文件并验证您是否拥有该网站。
如果找不到，请打开`https://xxxxx/googleexxxxxx.html`验证文件是否正确上传 


## 生成网站地图  

- **网页版:** [访问 XML Sitemaps](http://www.xml-sitemaps.com/)  
  
- **客户端:** [下载 Sitemap X](http://cn.sitemapx.com/)  

下载生成的地图文件 `sitemap.xml` 并上传至网站目录 `public/sitemap.xml`。
如果找不到，请打开`https://xxxxx/sitemap.xml`验证文件是否正确上传 


## robots.txt 文件  
写入txt文件内容
``` 
User-agent: *  
Disallow:  
Sitemap: https://sereinna.github.io/sitemap.xml
```
生成的roboes文件 `robots.txt` 并上传至网站根目录。


## 最后
通过以下链接手动添加新的站点地图：[添加新的站点地图](https://search.google.com/search-console/sitemaps)。  