---
title: markdown
published: 2024-09-16
description: "如何修改代码信息"
image: "./2.jpg"
tags: [linux,指南]
category: 指南
draft: false
---


# 如何在 Markdown 中展示代码  

Markdown以一种简单的方式来格式化代码。本文将介绍如何在 Markdown 中展示不同语言的代码。  

使用三个反引号（```）来创建代码块，并指定语言以启用语法高亮。 

### Bash 脚本  

```bash  
# Bash 代码示例  
echo "Hello, World!"  
```  

### Python 代码  

```python  
# Python 代码示例  
def greet():  
    print("Hello, World!")  

greet()  
```  

### JavaScript 代码  

```javascript  
// JavaScript 代码示例  
console.log("Hello, World!");  
```  

### YAML 配置  

```yaml  
# YAML 配置示例  
site: "https://your.domain"  
database:  
  host: "localhost"  
  port: 5432  
```  

### Diff 格式  

```diff  
- site: "https://fuwari.vercel.app/"  
+ site: "https://your.domain"  
```  

对于短代码片段，可以使用单个反引号（`）来展示内联代码。  

**示例**：  

要打印消息到控制台，你可以使用 `console.log("Hello, World!");`。  
