---
title: Fuwari的简单指南
published: 2024-09-16
description: "如何使用此博客模板。"
image: "./tiedao.jpeg"
tags: [Fuwari]
category: 示例
draft: false
---

> 封面图片来源: [来源](<https://image.civitai.com/xG1nkqKTMzGDvpLrqFT7WA/208fc754-890d-4adb-9753-2c963332675d/width=2048/01651-1456859105-(colour_1.5),girl,_Blue,yellow,green,cyan,purple,red,pink,_best,8k,UHD,masterpiece,male%20focus,%201boy,gloves,%20ponytail,%20long%20hair,.jpeg>)

这个博客模板是用 [Astro](https://astro.build/)。对于本指南中未提及的内容，您可以在 [Astro Docs](https://docs.astro.build/)。中找到答案

## Front-matter of Posts

```yaml
---
title: 我的第一篇博客文章
published: 2023-09-016
description: 这是我新Astro博客的第一篇文章。
image: ./cover.jpg
tags: [blog, write]
category: blog撰写
draft: false
---
```

| 属性          | 描述                                                                                                                                                               |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `title`       | 文章的标题。                                                                                                                                                       |
| `published`   | 文章发布的日期。                                                                                                                                                   |
| `description` | 文章的简短描述。在索引页上显示。                                                                                                                                   |
| `image`       | 文章的封面图片路径。<br/>1. 以 `http://` 或 `https://`开头：使用网络图片<br/>2. 以 `/`开头：用于 `public` 目录中的图片<br/>3. 没有前缀：相对于 Markdown 文件的位置 |
| `tags`        | 文章的标签。                                                                                                                                                       |
| `category`    | 文章的分类。                                                                                                                                                       |
| `draft`       | 如果这篇文章还是草稿，它将不会显示。                                                                                                                               |

## 文章文件应放置的位置

您的文章文件应该放置在 `src/content/posts/` 目录中。您还可以创建子目录来更好地组织您的文章和资产。

```
src/content/posts/
├── post-1.md
└── post-2/
    ├── cover.png
    └── index.md
```
