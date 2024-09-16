---
title: markdown
published: 2024-09-16
description: "如何修改代码信息"

tags: [Fuwari]
category: 指南
draft: false
---


# 如何在 Markdown 中展示代码  

在撰写技术博客或文档时，清晰地展示代码是非常重要的。Markdown 提供了一种简单的方式来格式化代码，使其更易读。本文将介绍如何在 Markdown 中展示不同语言的代码。  

使用三个反引号（```）来创建代码块，并指定语言以启用语法高亮。以下是一些常见的示例：  

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

通过使用适当的语言标记和格式，你可以确保代码在 Markdown 渲染器中正确显示，并且具有语法高亮功能。这不仅提高了代码的可读性，也帮助读者更容易地理解代码内容。