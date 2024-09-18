---
title: ollama
published: 2024-09-18
description: "如何部署ollama及其webui"
image: "https://www.pixiv.net/artworks/122436164"
tags: [win,ollama,model]
category: "指南"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 下载 `ollama`\
相关文章:[ollama](https://github.com/ollama/ollama?tab=readme-ov-file)

:::

##  下载并运行一个模型
```python 
ollama run xxx
```

##  更改win的下载路径
```python 
setx OLLAMA_MODELS "E:\ollama" /M
```

##  ollama命令行设置

| 属性          | 描述                                                                                                                                                               |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |

| `serve`    | 启动 ollama 服务。                                           |  
| `create`   | 使用 Modelfile 创建一个模型。这可能与定义模型架构或设置有关。|  
| `show`     | 显示特定模型的信息。                                         |  
| `run`      | 执行或运行指定的模型。                                       |  
| `pull`     | 从注册表中下载模型。                                         |  
| `push`     | 将模型上传到注册表。                                         |  
| `list`     | 列出所有模型。                                               |  
| `ps`       | 列出正在运行的模型。                                         |  
| `cp`       | 复制一个模型。                                               |  
| `rm`       | 删除一个模型。                                               |  
| `help`     | 获取关于任何命令的帮助。                                     |


##  ollama下载的模型与huggingface的模型的区别
通常情况下，Qwen模型的表示方法为Qwen1.5-4B-Chat。在Ollama中，Qwen指代的是与Hugging Face上的Qwen1_5-4B-Chat-q4_0.gguf模型相对应的版本，这是一个经过4位量化处理的模型

##  webui的使用