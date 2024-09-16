---
title: Markdown扩展功能
published: 2024-05-01
description: '了解更多关于Fuwari中的Markdown功能'
image: ''
tags: [Markdown, Fuwari]
category: '示例'
draft: true 
---

## GitHub仓库卡片
您可以添加动态卡片，链接到GitHub仓库。页面加载时，仓库信息将从GitHub API中提取。

::github{repo="Fabrizz/MMM-OnSpotify"}

使用代码 `::github{repo="<所有者>/<仓库>"}` 创建一个GitHub仓库卡片。

```markdown
::github{repo="saicaca/fuwari"}
```

## 警告

以下类型的警告得到支持：`note` `tip` `important` `warning` `caution`

:::note
突出显示用户在浏览时应考虑的信息。
:::

:::tip
提供帮助用户更成功的可选信息。
:::

:::important
对于用户成功至关重要的关键信息。
:::

:::warning
由于潜在风险，需要立即引起用户注意的关键内容。
:::

:::caution
行动的负面潜在后果。
:::

```markdown
:::note
突出显示用户在浏览时应考虑的信息。
:::

:::tip
提供帮助用户更成功的可选信息。
:::
```

警告的标题可以自定义。

:::note[我的自定义标题]
这是一个具有自定义标题的笔记。
:::

```markdown
:::note[我的自定义标题]
这是一个具有自定义标题的笔记。
:::
```

> [!TIP]
> [The GitHub syntax](https://github.com/orgs/community/discussions/16925) 也得到了支持。

```
> [!NOTE]
> GitHub语法也得到了支持。

> [!TIP]
> GitHub语法也得到了支持。
```