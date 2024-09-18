---
title: ollama
published: 2024-09-18
description: "如何部署ollama及其webui"
image: "./122436164_p0.png"
tags: [win,ollama,model]
category: "指南"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 下载 `ollama`\
相关文章:  [ollama](https://github.com/ollama/ollama?tab=readme-ov-file)

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
| `serve`     | 启动 ollama 服务。                                                                      |
| `create`    | 使用 Modelfile 创建一个模型。这可能与定义模型架构或设置有关。                            |
| `show`      | 显示特定模型的信息。                                                                    |
| `run`       | 执行或运行指定的模型。                                                                  |
| `pull`      | 从注册表中下载模型。                                                                    |
| `push`      | 将模型上传到注册表。                                                                    |
| `list`      | 列出所有模型。                                                                          |
| `ps`        | 列出正在运行的模型。                                                                    |
| `cp`        | 复制一个模型。                                                                          |
| `rm`        | 删除一个模型。                                                                          |
| `help`      | 获取关于任何命令的帮助。                                                                |


##  ollama下载的模型与huggingface的模型的区别
通常情况下，Qwen模型的表示方法为Qwen1.5-4B-Chat
在Ollama中，Qwen指代的是与Hugging Face上的Qwen1_5-4B-Chat-q4_0.gguf模型相对应的版本，这是一个经过4位量化处理的模型



##  ollama模型列表
| Model                | Parameters | Size  | Download                       |
|----------------------|------------|-------|---------------------------------|
| Llama 3.1            | 8B         | 4.7GB | `ollama run llama3.1`           |
| Llama 3.1            | 70B        | 40GB  | `ollama run llama3.1:70b`       |
| Llama 3.1            | 405B       | 231GB | `ollama run llama3.1:405b`      |
| Phi 3 Mini           | 3.8B       | 2.3GB | `ollama run phi3`               |
| Phi 3 Medium         | 14B        | 7.9GB | `ollama run phi3:medium`        |
| Gemma 2              | 2B         | 1.6GB | `ollama run gemma2:2b`          |
| Gemma 2              | 9B         | 5.5GB | `ollama run gemma2`             |
| Gemma 2              | 27B        | 16GB  | `ollama run gemma2:27b`         |
| Mistral              | 7B         | 4.1GB | `ollama run mistral`            |
| Moondream 2          | 1.4B       | 829MB | `ollama run moondream`          |
| Neural Chat          | 7B         | 4.1GB | `ollama run neural-chat`        |
| Starling             | 7B         | 4.1GB | `ollama run starling-lm`        |
| Code Llama           | 7B         | 3.8GB | `ollama run codellama`          |
| Llama 2 Uncensored   | 7B         | 3.8GB | `ollama run llama2-uncensored`  |
| LLaVA                | 7B         | 4.5GB | `ollama run llava`              |
| Solar                | 10.7B      | 6.1GB | `ollama run solar`              |


##  webui的使用(之后再说)

### way_1安装docker
[docker](https://desktop.docker.com/win/stable/amd64/Docker%20Desktop%20Installer.exe)


### way_2python安装

对于 Windows：

```bash
git clone https://github.com/open-webui/open-webui.git
cd open-webui

copy .env.example .env

npm install
npm run build

cd .\backend
```

使用 Conda 作为开发环境进行安装：

```bash
# 创建并激活 Conda 环境
conda create --name open-webui-env python=3.11
conda activate open-webui-env
```

执行以下命令安装依赖：

```bash
pip install -r requirements.txt -U
```

运行启动脚本：

```bash
 .\start_windows.bat
```

在 [http://localhost:8080/](http://localhost:8080/) 上启动并运行 Open WebUI 

