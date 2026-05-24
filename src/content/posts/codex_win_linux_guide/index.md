---
title: Codex 在 Windows 和 Linux 上的调用指南
published: 2026-05-24
description: "一篇写清楚 Codex CLI、VS Code 扩展、config.toml、Windows 原生与 Linux/WSL 使用方式的入门指南。"
image: "./145065828_p0.jpg"
tags: [Codex, Windows, Linux, VSCode, 指南]
category: "AI工具"
draft: false
lang: ""
---

![插图1](./145065828_p0.jpg)

## 前言

:::tip
如果你主要在 Linux 上工作，当前最稳的入口是 **Codex CLI** 和 **VS Code 侧边栏扩展**。  
如果你在 Windows 上工作，优先推荐 **PowerShell 原生运行**；只有在你依赖 Linux 工具链时，再切换到 **WSL2**。
:::

根据 OpenAI 官方 Quickstart，Codex 桌面 App 页面目前提供 **Windows 下载**，Linux 仍显示 **Get notified for Linux**。这意味着在 Linux 上，现阶段更适合直接使用 CLI 和 IDE 扩展，而不是等待原生桌面版。

## 一、Codex 目前有哪些入口

常见入口有 3 个：

1. `codex` 命令行
2. VS Code/Cursor/Windsurf 侧边栏扩展
3. Windows 桌面 App

如果你的目标是“能写代码、能读仓库、能跑命令、能改文件”，最实用的还是前两种：

- **CLI** 适合终端党、远程机器、Linux 服务器、WSL
- **VS Code 扩展** 适合边看代码边提问、边选中文件边让 Codex 修改

## 二、Windows 和 Linux 的安装方式

### 1. 命令行安装 Codex

OpenAI 官方 CLI 文档当前给出的安装命令是：

```bash
npm i -g @openai/codex
```

首次运行：

```bash
codex
```

第一次启动时，Codex 会提示你登录，可以使用：

- ChatGPT 账号登录
- OpenAI API Key 登录

### 2. Windows 安装示例

先确认已经安装 `Node.js` 和 `npm`，然后在 PowerShell 中执行：

```powershell
npm i -g @openai/codex
codex
```

常用补充命令：

```powershell
codex login
codex doctor
codex update
```

如果你希望升级到最新版：

```powershell
npm i -g @openai/codex@latest
```

### 3. Linux 安装示例

在 Ubuntu、Debian、Arch、Fedora 一类 Linux 环境中，本质也是一样的：

```bash
npm i -g @openai/codex
codex
```

常用补充命令：

```bash
codex login
codex doctor
codex update
```

### 4. Linux 不走 npm 的两种办法

如果你不想用 `npm`，Linux 里还可以走两条路：

#### 方案 A：直接 `git clone` 源码构建

这是 OpenAI 官方 `openai/codex` 仓库 `docs/install.md` 明确写出的方式，适合：

- 不想依赖全局 `npm`
- 想自己编译最新版
- 想顺便阅读或修改 Codex 源码

```bash
git clone https://github.com/openai/codex.git
cd codex/codex-rs

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

rustup component add rustfmt
rustup component add clippy
cargo install --locked just
cargo install --locked cargo-nextest

cargo build
cargo run --bin codex -- "explain this codebase to me"
```

如果你希望把它装到自己的 Cargo bin 目录，通常也可以进一步执行：

```bash
cargo install --path .
```

安装后确认：

```bash
~/.cargo/bin/codex --version
```

#### 方案 B：下载 GitHub Release 压缩包后手动解压

如果你不想编译，也可以直接从官方 GitHub Releases 下载 Linux 包，再手动解压。

你这个流程更准确，核心就是：

- 下载官方发布的 `tar.gz`
- 本地解压
- 把解压出的可执行文件移动到 `~/.local/bin`
- 加入 `PATH` 后直接使用

推荐写法如下：

```bash
curl -L -o codex.tar.gz \
  https://github.com/openai/codex/releases/latest/download/codex-x86_64-unknown-linux-musl.tar.gz

tar -xzf codex.tar.gz

mkdir -p ~/.local/bin
mv codex-x86_64-unknown-linux-musl ~/.local/bin/codex
chmod +x ~/.local/bin/codex

echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

codex --version
```

这条路线更像“下载二进制并解压即用”，适合：

- 机器上不想装 Node.js
- 也不想装 Rust 工具链
- 只想快速拿到一个可执行文件

要注意两点：

- 上面示例以 `x86_64` Linux 为例，不同架构要换成对应发布包
- GitHub Release 资产命名会随平台变化，下载前最好先看一眼 Releases 页面

### 5. Windows 到底选 PowerShell 还是 WSL2

OpenAI 官方 Windows 文档明确建议：

- **默认优先使用 Windows 原生运行**
- **需要 Linux 原生工具链时再用 WSL2**

也就是说：

- 你的项目在 `C:\`、工具主要是 PowerShell / Windows Git / Windows Node.js，就直接原生跑
- 你的项目在 WSL2 里，依赖 `bash`、`apt`、`python venv`、`docker`、`make` 等 Linux 工作流，就把 Codex 切到 WSL

## 三、VS Code 安装 Codex

### 1. 图形界面安装

打开 VS Code：

1. 进入扩展市场
2. 搜索 `Codex`
3. 安装 OpenAI 官方扩展 `Codex - OpenAI's coding agent`

安装完成后，如果没有立刻看到侧边栏，官方文档建议 **重启 VS Code**。

### 2. 命令行安装扩展

如果你的 VS Code 命令行可用，也可以直接安装：

```bash
code --install-extension OpenAI.chatgpt
```

这个扩展当前在 Marketplace 的条目 ID 是 `OpenAI.chatgpt`。

### 3. 扩展装好后你会看到什么

根据官方 IDE 文档：

- Codex 会出现在编辑器侧边栏
- 在 VS Code 中默认会打开在**右侧边栏**
- 它和 CLI 使用的是**同一个 agent**
- 它和 CLI 共享**同一套配置**

## 四、`config.toml` 放哪里，CLI 和侧边栏共不共用

答案是：**共用同一份配置**。

官方 Config Basics 文档说明：

- 用户级配置文件：`~/.codex/config.toml`
- 项目级覆盖配置：`.codex/config.toml`
- CLI 和 IDE Extension 共享同样的配置层

### 1. Windows 路径

```powershell
$env:USERPROFILE\.codex\config.toml
```

通常就是：

```text
C:\Users\你的用户名\.codex\config.toml
```

### 2. Linux 路径

```bash
~/.codex/config.toml
```

### 3. 在 VS Code 侧边栏里直接打开配置

官方文档给出的路径是：

1. 点击 Codex 侧边栏右上角齿轮
2. 选择 `Codex Settings`
3. 选择 `Open config.toml`

### 4. 配置优先级

从高到低大致是：

1. CLI 启动参数和 `--config`
2. `--profile <name>`
3. 项目内 `.codex/config.toml`
4. 用户级 `~/.codex/config.toml`
5. 系统配置
6. 内置默认值

## 五、推荐直接可用的 `config.toml`

下面这份更适合作为个人入门配置。Windows 用户可以保留 `[windows]`，Linux 用户删掉这一段即可。

```toml
model = "gpt-5.5"
model_provider = "openai"
approval_policy = "on-request"
sandbox_mode = "workspace-write"
web_search = "cached"
model_reasoning_effort = "medium"
model_verbosity = "medium"
# personality = "pragmatic"

[profiles.safe]
approval_policy = "untrusted"
sandbox_mode = "read-only"

[profiles.full]
approval_policy = "never"
sandbox_mode = "danger-full-access"

[windows]
sandbox = "elevated"
```

### 这几个字段最常用

| 字段 | 作用 | 建议 |
| --- | --- | --- |
| `model` | 默认模型 | 先用官方 sample 里的 `gpt-5.5` |
| `approval_policy` | 什么时候需要你批准 | 日常推荐 `on-request` |
| `sandbox_mode` | 文件和命令权限范围 | 日常推荐 `workspace-write` |
| `web_search` | 是否启用搜索 | 一般用 `cached` |
| `model_reasoning_effort` | 思考强度 | 日常 `medium` 即可 |
| `model_verbosity` | 输出啰嗦程度 | `medium` 很稳 |
| `[windows].sandbox` | Windows 原生沙箱模式 | 官方推荐 `elevated` |

### 适合新手的两个原则

- 不要一上来就用 `danger-full-access`
- 不要把 `approval_policy = "never"` 当默认配置

更稳的方式是：

- 平时用 `on-request + workspace-write`
- 确认环境安全时，再临时切到更高权限

## 六、命令行如何使用 Codex

### 1. 直接进入交互模式

```bash
codex
```

进入后你可以像聊天一样直接提需求，例如：

```text
帮我阅读这个仓库，并总结 src 目录结构
```

或者：

```text
修复当前测试失败，但不要改 public API
```

### 2. 直接带一句 prompt 启动

```bash
codex "解释这个仓库的构建流程"
```

### 3. 非交互执行：`codex exec`

适合脚本、批处理、CI：

```bash
codex exec "summarize the repo structure"
```

允许在工作区内改文件：

```bash
codex exec --sandbox workspace-write "fix the failing test"
```

输出 JSONL，适合脚本消费：

```bash
codex exec --json "summarize the repo structure"
```

把最终回答单独写入文件：

```bash
codex exec "Extract project metadata" -o ./codex-output.txt
```

### 4. 常见 CLI 子命令

```bash
codex review
codex doctor
codex update
codex mcp list
codex mcp --help
```

其中：

- `codex review`：让 Codex 做一次代码审查
- `codex doctor`：检查安装、登录、配置、运行环境
- `codex update`：升级 Codex
- `codex mcp list`：查看 MCP 服务

### 5. 推荐记住的 CLI 斜杠命令

交互界面里最常用的是这些：

```text
/status
/model
/permissions
/diff
/review
/init
/mcp
```

简要理解如下：

- `/status`：看当前模型、权限、上下文占用
- `/model`：切换模型与推理强度
- `/permissions`：快速切换权限模式
- `/diff`：查看本次改动的 Git diff
- `/review`：让 Codex 复审当前工作区
- `/init`：生成 `AGENTS.md`
- `/mcp`：查看当前 MCP 工具是否可用

## 七、侧边栏如何使用 Codex

### 1. 最基础的工作流

打开侧边栏后，一般按这个顺序用就够了：

1. 选择当前项目目录
2. 确认模式是 `Local`
3. 在输入框里直接描述任务
4. 需要引用文件时，用 `@文件名`
5. 看 diff、看回答、决定是否继续追问

### 2. 如何把文件上下文喂给 Codex

官方 IDE Features 文档支持两种很实用的方式：

- 直接在 prompt 里写 `@example.tsx`
- 选中代码后，把选区加入当前线程

示例：

```text
参考 @src/pages/index.tsx 的风格，新增一个 Resources 页面
```

### 3. 侧边栏里几个最常用的开关

#### 模型切换

输入框下方的模型切换器可以改：

- 模型
- reasoning effort（如 `low`、`medium`、`high`）

#### 模式切换

常见有三种理解：

- `Chat`：只讨论、只规划，不急着动文件
- `Agent`：允许读代码、改文件、运行命令
- `Agent (Full Access)`：高权限模式，谨慎使用

### 4. VS Code 命令面板里可直接调用的命令

按 `Ctrl+Shift+P` 后搜索 `Codex`，常用命令包括：

- `chatgpt.openSidebar`
- `chatgpt.newChat`
- `chatgpt.addToThread`
- `chatgpt.addFileToThread`

### 5. 侧边栏里的斜杠命令

官方 IDE 文档列出的常用命令包括：

```text
/auto-context
/local
/cloud
/review
/status
/goal
```

它们分别适合：

- `/auto-context`：自动带上最近文件和 IDE 上下文
- `/local`：切回本地执行
- `/cloud`：切到云端任务
- `/review`：开启代码审查模式
- `/status`：查看线程、上下文、限制信息
- `/goal`：为当前任务设一个持续目标

## 八、Windows 专属设置：侧边栏走原生还是 WSL

如果你在 Windows 上使用 VS Code 侧边栏，官方 Settings 文档里有一个非常关键的设置：

```text
chatgpt.runCodexInWindowsSubsystemForLinux
```

含义是：

- `false`：Codex 在 Windows 原生环境运行
- `true`：Codex 在 WSL 中运行

什么时候改成 `true`：

- 你的仓库在 WSL 文件系统里
- 你的依赖管理、Python 环境、Docker、构建脚本都在 Linux 下
- 你明确需要 Linux 原生命令

什么时候保持默认原生 Windows：

- 你的项目在 Windows 目录中
- 你主要用 PowerShell、Windows Git、Windows Node.js
- 你更在意原生速度和更直接的 Windows 沙箱

## 九、MCP 怎么接进来

如果你想让 Codex 连接额外工具，官方 CLI 方式是：

```bash
codex mcp add context7 -- npx -y @upstash/context7-mcp
```

更通用的写法：

```bash
codex mcp add <server-name> --env VAR1=VALUE1 -- <stdio-server-command>
```

查看是否生效：

```bash
codex mcp list
```

在 CLI 交互界面里也可以：

```text
/mcp
```

## 十、一个简单的上手顺序

如果你是第一次装 Codex，我建议直接按下面顺序来：

1. 安装 Node.js / npm
2. `npm i -g @openai/codex`
3. 运行 `codex`
4. 登录 ChatGPT 账号或 API Key
5. 写好 `~/.codex/config.toml`
6. 在 VS Code 安装 `OpenAI.chatgpt`
7. 进入侧边栏，用 `@文件名 + 自然语言` 开始第一次任务

## 十一、我自己更推荐的用法

### 如果你是 Windows 用户

- 默认：**PowerShell 原生 Codex**
- 代码编辑：**VS Code 侧边栏**
- Linux 工具链项目：**切 WSL2**

### 如果你是 Linux 用户

- 默认：**CLI + VS Code 扩展**
- 项目级规则：加 `.codex/config.toml`
- 自动化：优先 `codex exec`

## 参考资料

- OpenAI Codex Quickstart: https://developers.openai.com/codex/quickstart
- OpenAI Codex CLI: https://developers.openai.com/codex/cli
- OpenAI Codex IDE: https://developers.openai.com/codex/ide
- OpenAI Codex Config Basics: https://developers.openai.com/codex/config-basic
- OpenAI Codex Config Reference: https://developers.openai.com/codex/config-reference
- OpenAI Codex Windows: https://developers.openai.com/codex/windows
- OpenAI Codex Authentication: https://developers.openai.com/codex/auth
- OpenAI Codex MCP: https://developers.openai.com/codex/mcp
- VS Code Marketplace: https://marketplace.visualstudio.com/items?itemName=OpenAI.chatgpt
- OpenAI Codex install.md: https://github.com/openai/codex/blob/main/docs/install.md
- OpenAI Codex Releases: https://github.com/openai/codex/releases

![插图2](./145089401_p0.jpg)
