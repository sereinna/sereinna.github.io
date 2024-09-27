---
title: Claude修饰分子的论文
published: 2024-09-27
description: "作为分子设计引擎的大型语言模型Large Language Models asMolecular Design Engines"
image: "./101325591_p0_master1200.jpg"
tags: [gpt,smiles,文献阅读]
category: "论文"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 文献阅读器 `zetoro`，
[相关文献](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01396)
[相关公众号](https://mp.weixin.qq.com/s/axtTU_qz51H6hcVx_IRGcA)



## 摘要
预训练的大型语言模型（LLM）已成为分子设计的潜在工具，因为它们似乎能够根据自然语言提示提供的简单指令来创建和修改分子。这项工作展示了 Claude3OpusLLM 可以根据提示读取、写入和修改分子，其中有效且独特的分子高达 97%

当被要求使用简单的自然语言提示操纵分子的电子结构时，该模型能够执行引导分子生成

## 介绍
### 研究背景与目的:
小分子设计对从药物发现到能源存储等多个领域至关重要。
传统的生成式机器学习方法存在训练复杂、生成分子有效性低等问题。
研究旨在探索预训练LLM在分子设计中的应用潜力。

### 研究方法:
使用Claude 3 Opus模型进行分子生成和修改实验。
采用SMILES字符串表示分子结构。
通过不同的自然语言提示词来引导模型生成分子。
使用摩根指纹和PCA构建分子的低维潜在空间表示。

### 主要发现:
LLM能够根据提示词读取、写入和修改分子,生成的分子97%是有效且独特的。
模型能够根据提示词控制生成分子与原始分子的相似度。
模型成功实现了引导式分子生成,如添加电子给体或吸收体基团。
通过分析潜在空间中的位移,揭示了不同提示词对分子生成的影响。

### 结论与展望:
LLM展现出作为强大且多功能的分子设计引擎的潜力。
相比传统方法,LLM具有更高的灵活性和泛化能力。
未来研究可以探索自动提示工程等方法,进一步提高LLM在分子设计中的应用效果。

### 研究意义:
为使用简单自然语言指令进行分子设计提供了新思路。
有望加速具有特定性质的新型分子的设计过程。
总的来说,这项研究展示了大型语言模型在分子设计领域的应用前景,为未来的药物发现和材料设计等领域提供了新的研究方向。


## 处理步骤

### 使用的模型
- 文章使用了Anthropic的Claude 3 Opus模型进行实验。

### 与模型交互的方式
- 使用Anthropic Python SDK进行交互。
- 请求包含任务指令（提示词）被发送到Anthropic服务器上的预训练模型进行处理。

### 模型参数设置
- temperature设置为0，以获得更确定性和集中的输出。
- max_tokens设置为1024，限制生成输出的长度。

### 基础系统提示词
对于每个查询，都提供以下系统提示词：
"You are a chemoinformatics expert that can generate new molecules. Please provide only the Python formatted list of SMILES strings, like [SMILES1, SMILES2, SMILES3] without any additional explanations or text."

### 任务提示词模板
"Given the molecule with SMILES representation 'smiles', generate n molecules that are prompt_detail. Respond with just the SMILES strings as elements of a Python list."
其中：
- smiles被替换为父分子的SMILES表示
- n被替换为目标候选分子数量（10）
- prompt_detail被替换为特定任务描述

### 具体提示词示例
文章中使用了多种提示词，包括：
- 基础提示词（A-H）：用于生成相似或完全不同的分子
- 引导生成提示词（I-P）：用于添加电子给体或吸收体基团
- 控制生成提示词（Q-S）：用于控制生成分子与父分子的相似度

### 模型输出
- 模型返回一个Python列表格式的SMILES字符串，代表生成的分子。

### 输出处理
- 使用RDKit验证返回的SMILES字符串的有效性。
- 将有效分子转换为规范形式。
- 过滤掉重复项和未修改的父分子。

### 评估指标
- Tanimoto相似度：评估生成分子与父分子的相似程度。
- 有效率：评估生成的有效且独特分子的比例。
- 化学多样性：评估生成分子集合的异质性。

### 模型性能
- 97%的输出是有效且独特的分子。
- API调用的中位响应时间为10.4秒（生成10个分子）。

### 可视化
- 研究者开发了一个基于Python的查看器，用于显示输入分子和每个提示词生成的输出分子。



## prompt
### 基础提示词 (A-H):

| 属性 | 描述 | 中文翻译 |
|------|------|----------|
| A | similar molecules by changing one or two atoms or bonds to produce closely related structures | 通过改变一两个原子或键来产生密切相关的结构，生成相似分子 |
| B | similar molecules by tweaking only the side chains | 仅通过调整侧链来生成相似分子 |
| C | similar molecules with minimal structural changes to find similar but new candidates | 通过最小的结构变化来寻找相似但新的候选分子 |
| D | similar molecules with slight variations on functional groups while maintaining the backbone structure | 在保持骨架结构的同时，对官能团进行轻微变化，生成相似分子 |
| E | completely different molecules by changing multiple atoms or bonds | 通过改变多个原子或键来生成完全不同的分子 |
| F | completely different molecules by significantly altering the core structure and introducing completely new functional groups | 通过显著改变核心结构并引入全新的官能团来生成完全不同的分子 |
| G | completely different molecules that significantly vary in size and functional groups | 生成在大小和官能团方面有显著差异的完全不同分子 |
| H | completely different molecules with significant structural changes to find new candidates | 通过显著的结构变化来寻找新的候选分子，生成完全不同的分子 |

### 引导生成提示词 (I-P):

| 属性 | 描述 | 中文翻译 |
|------|------|----------|
| I | Similar molecules by changing one or two atoms or bonds to produce closely related structures focusing on incorporating electron donating groups (EDGs) to find new candidates | 通过改变一两个原子或键来产生密切相关的结构，重点是引入电子给体基团(EDG)以寻找新的候选分子 |
| J | Similar molecules by tweaking only the side chains to produce closely related structures focusing on incorporating electron donating groups (EDGs) to find new candidates | 仅通过调整侧链来产生密切相关的结构，重点是引入电子给体基团(EDG)以寻找新的候选分子 |
| K | Similar molecules with minimal structural changes to produce closely related structures focusing on incorporating electron donating groups (EDGs) to find new candidates | 通过最小的结构变化来产生密切相关的结构，重点是引入电子给体基团(EDG)以寻找新的候选分子 |
| L | Similar molecules with slight variations on functional groups while maintaining the backbone structure to produce closely related structures focusing on incorporating electron donating groups (EDGs) to find new candidates | 在保持骨架结构的同时，对官能团进行轻微变化来产生密切相关的结构，重点是引入电子给体基团(EDG)以寻找新的候选分子 |
| M | Similar molecules by changing one or two atoms or bonds to produce closely related structures focusing on incorporating electron withdrawing groups (EWGs) to find new candidates | 通过改变一两个原子或键来产生密切相关的结构，重点是引入电子吸引基团(EWG)以寻找新的候选分子 |
| N | Similar molecules by tweaking only the side chains to produce closely related structures focusing on incorporating electron withdrawing groups (EWGs) to find new candidates | 仅通过调整侧链来产生密切相关的结构，重点是引入电子吸引基团(EWG)以寻找新的候选分子 |
| O | Similar molecules with minimal structural changes to produce closely related structures focusing on incorporating electron withdrawing groups (EWGs) to find new candidates | 通过最小的结构变化来产生密切相关的结构，重点是引入电子吸引基团(EWG)以寻找新的候选分子 |
| P | Similar molecules with slight variations on functional groups while maintaining the backbone structure to produce closely related structures focusing on incorporating electron withdrawing groups (EWGs) to find new candidates | 在保持骨架结构的同时，对官能团进行轻微变化来产生密切相关的结构，重点是引入电子吸引基团(EWG)以寻找新的候选分子 |

### 控制生成提示词 (Q-S):

| 属性 | 描述 | 中文翻译 |
|------|------|----------|
| Q | Barely similar (very low Tanimoto similarity) molecules compared to the given parent molecule by altering atoms, bonds, functional groups, or making other changes to find new candidates | 与给定的父分子相比，通过改变原子、键、官能团或进行其他变化来寻找几乎不相似（非常低的Tanimoto相似度）的新候选分子 |
| R | Marginally similar (low Tanimoto similarity) molecules compared to the given parent molecule by altering atoms, bonds, functional groups, or making other changes to find new candidates | 与给定的父分子相比，通过改变原子、键、官能团或进行其他变化来寻找略微相似（低Tanimoto相似度）的新候选分子 |
| S | Moderately similar (moderate Tanimoto similarity) molecules compared to the given parent molecule by altering atoms, bonds, functional groups, or making other changes to find new candidates | 与给定的父分子相比，通过改变原子、键、官能团或进行其他变化来寻找中等相似（中等Tanimoto相似度）的新候选分子 |
这些提示词都是在基础系统提示词的基础上使用的。完整的提示词结构如下：


### 基础系统提示词
系统提示词：
"You are a chemoinformatics expert that can generate new molecules. Please provide only the Python formatted list of SMILES strings, like [SMILES1, SMILES2, SMILES3] without any additional explanations or text."
任务提示词模板：
"Given the molecule with SMILES representation 'smiles', generate n molecules that are prompt_detail. Respond with just the SMILES strings as elements of a Python list."
其中，'smiles'会被替换为特定的父分子SMILES表示，n会被替换为10（要生成的分子数量），prompt_detail会被替换为上述A-S中的某一个具体提示词内容。

这种结构允许研究者系统地测试不同类型的分子修改指令，从而评估Claude模型在分子设计任务中的表现和能力。