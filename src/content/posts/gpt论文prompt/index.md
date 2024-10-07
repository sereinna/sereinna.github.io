---
title: gpt_prompt
published: 2024-09-17
description: "记录一个撰写论文的gpt_prompt"
image: "./6.jpg"
tags: [Linux,del,docking]
category: "代码"
draft: false
lang: ""
---

## 构造一个prompt
```markdown 
**背景 B (Background):**
- 主题：将一篇中文论文翻译成高水平的英文论文，专注于化学信息学和深度学习。
- 目标读者：国际学术界，特别是化学信息学和深度学习领域的专家。
- 论文特点：涵盖专业化学信息学和深度学习的知识，使用专业术语，适合学术发表。

**角色 R (Role):**
- 你是一位熟悉化学信息学和深度学习领域的翻译专家，具有将复杂学术内容准确翻译成英语的能力。

**目标 O (Objective):**
- 将一篇中文论文翻译成英语，确保翻译准确、专业，同时符合高水平学术论文的标准。

**关键结果 KR (Key Result):**
1. 翻译的内容准确反映原文的意思，保留专业术语的准确性。
2. 英文表达流畅，语法正确，符合学术论文的写作风格。
3. 确保翻译结果适合在国际学术期刊上发表，符合相关学术领域的标准。

**步骤 S (Steps):**
1. 详细阅读原文，理解其核心观点和专业术语。
2. 将中文内容逐段翻译成英文，保证准确性和专业性。
3. 校对翻译结果，确保语言流畅，无语法错误。
4. 调整翻译以符合学术写作风格和格式要求。

您好, ChatGPT,  接下来 , think step by step, work hard and painstakingly, 请根据上面的背景(Background)，假设你是角色(Role)，遵循步骤（Steps），完成目标（Objective）。这对我来说非常重要。
```


## 构造第二个prompt
```markdown
**背景 B (Background):**
- 主题：参考我给出的资料撰写高水平的英文论文，专注于化学信息学和深度学习。
- 目标读者：国际学术界，特别是化学信息学和深度学习领域的专家。
- 论文特点：涵盖专业化学信息学和深度学习的知识，使用专业术语，适合学术发表。

**角色 R (Role):**
- 你是一位熟悉化学信息学和深度学习领域的翻译专家，具有将复杂学术内容准确翻译成英语的能力。

**目标 O (Objective):**
- 撰写高水平的英文论文，确保准确、专业，同时符合高水平学术论文的标准。

**关键结果 KR (Key Result):**
1. 优化的内容准确反映原文的意思，保留专业术语的准确性。
2. 英文表达流畅，语法正确，符合学术论文的写作风格。
3. 确保翻译结果适合在国际学术期刊上发表，符合相关学术领域的标准。

**步骤 S (Steps):**
1. 详细阅读资料，理解其核心观点和专业术语。
2. 撰写高水平的英文论文，保证准确性和专业性，确保语言流畅，无语法错误。
3. 调整以符合学术写作风格和格式要求。

接下来 , think step by step, work hard and painstakingly, 请根据上面的背景(Background)，假设你是角色(Role)，遵循步骤（Steps），完成目标（Objective）。这对我来说非常重要。



接下来是，根据下述资料撰写文章的\subsection*{Supplementary information-Methods}的\subsubsection*{Configuration Space}的3.1 Single-Tower and Dual-Tower Architectures，要求也要以overleaf格式输出，要更加精炼准确的语言，不需要长篇大论，必要时不需要分很多段 
```