---
title: PyMOL蛋白质结构可视化完整指南
published: 2025-01-17
description: "详细的PyMOL操作指南，包含基础设置到高级可视化的完整教程"
image: "./125226669_p0.png"
tags: [pymol,图]
category: "指南"
draft: false
lang: ""
---


:::tip
准备工作 蛋白质 `pdb`，配体 `mol2`，画图软件 `pymol`
:::

# 前言

本指南提供了使用PyMOL进行蛋白质结构可视化的完整操作流程，从基础加载到高级渲染都有详细说明。所有代码块都可以直接复制使用，注意替换相关文件名和路径。

# 1. 基础操作

## 1.1 结构加载
结构加载是PyMOL使用的第一步，支持多种分子结构文件格式。

```bash
# 加载蛋白质结构文件
# 语法: load "文件路径", 对象名称
load "3p3h_protein.pdb", 3p3h_protein  # 加载PDB格式的蛋白质结构
load "docked_456.mol2", docked_456     # 加载MOL2格式的配体结构

# 创建结构副本
# 用于创建同一结构的不同表示方式
create 3p3h_protein2, 3p3h_protein     # 创建蛋白质结构的副本
```

## 1.2 基本显示模式
PyMOL提供多种分子显示模式，可以组合使用创建更好的可视化效果。

```bash
# 设置不同的显示模式
show cartoon, 3p3h_protein    # 显示蛋白质的卡通模式
show surface, 3p3h_protein2   # 显示蛋白质的表面
show sticks, docked_456       # 显示分子的棍状模型

# 隐藏特定显示模式
hide lines, all              # 隐藏所有对象的线框显示
hide everything, all         # 隐藏所有显示模式
```

# 2. 视觉参数设置

## 2.1 颜色调整
合理的颜色设置能够突出重要结构并提高可视化效果。

```bash
# 基础颜色设置
color white, 3p3h_protein     # 将蛋白质设为白色
color yellow, docked_456      # 将配体设为黄色

# 自定义颜色（RGB值范围0-1）
set_color custom_color1,[0.851, 0.867, 0.894] # 创建自定义颜色
color custom_color1, 5fl4_protein2              # 应用自定义颜色

# 特定结构着色
util.cba(36, "docked_1_1", _self=cmd)    # 使用色谱着色方案
util.cba(154, "obj02", _self=cmd)         # 对特定对象使用色谱
```

## 2.2 透明度和显示效果
调整透明度可以同时展示多个结构层次。

```bash
# 透明度设置（0完全透明-1完全不透明）
set transparency, 0.5, 3p3h_protein2           # 设置表面透明度
set cartoon_transparency, 0.7, 6w8i_pre        # 设置卡通图透明度

# 显示效果优化
set cartoon_fancy_helices, 1                   # 启用精美螺旋效果
set cartoon_gap_cutoff, 0                      # 设置卡通图断裂阈值
set stick_radius, 0.2 , xx                       # 设置棍状模型粗细
bg_color white                                 # 设置背景颜色为白色
```

## 2.3 标签和字体设置
清晰的标签对于结构说明至关重要。

```bash
# 字体设置
set label_font_id, 1    # 设置字体类型(0:位图 1:Sans 2:Serif 3:Courier)
set label_size, 20      # 设置字体大小

# 残基标签设置
label my_selection and name CA, "%s%s" % (resi,resn)  # 显示残基编号和名称
label site and name CA, "%s-%s" % (resi,single[resn]) # 显示简写格式(如"23-A")

# 单氨基酸代码定义
single={'VAL':'V','ILE':'I','LEU':'L','GLU':'E','GLN':'Q',
        'ASP':'D','ASN':'N','HIS':'H','TRP':'W','PHE':'F',
        'TYR':'Y','ARG':'R','LYS':'K','SER':'S','THR':'T',
        'MET':'M','ALA':'A','GLY':'G','PRO':'P','CYS':'C'}
```

# 3. 结构分析工具

## 3.1 残基选择
精确的残基选择是结构分析的基础。

```bash
# 基本选择语法
select residue_d88, chain D and resi 88                # 选择特定链的特定残基
select binding_site, (resi 263+267+271+128+132+198)    # 选择多个残基

# 基于距离的选择
select interacting_residues, (3p3h_protein within 4 of docked_456)  # 选择4Å范围内的残基

# 复杂选择条件
select selected_residues, (chain D and resi 88 and resn D) or (chain D and resi 140 and resn N)
```
## 3.2 结构分析
用于分析分子间相互作用和结构特征。
```bash
# 氢键分析
distance hbond, (resn MG), (resi 263+267+271+128+132+198), mode=2, cutoff=3.5

# 结构对齐
cealign 1,18                # 对齐蛋白结构
align mobile, target        # 对齐分子

# 相互作用分析
# 计算氢键（距离小于3.2Å）
distance hydrogen_bonds, (docked_456 and (name O+N)), (interacting_residues and (name O+N)), cutoff=3.2

# 计算盐桥（距离小于4.0Å）
distance salt_bridges, (docked_456 and (resn ARG+LYS)), (3p3h_protein and (resn ASP+GLU)), cutoff=4.0
```
# 4. 场景制作
## 4.1 远景图制作
展示整体结构的完整流程。
```bash
# 1. 基础设置
load "protein.pdb", protein
create protein_surface, protein

# 2. 显示设置
show cartoon, protein
show surface, protein_surface
show sticks, ligand

# 3. 视觉优化
color white, protein
color white, protein_surface
set transparency, 0.5, protein_surface
set cartoon_fancy_helices, 1
bg_color white

# 4. 渲染和保存
ray 1000, 1000
png "distant_view.png", dpi=300
```
## 4.2 近景图制作
展示具体相互作用的详细流程。
```bash
# 1. 选择和显示相互作用区域
select interacting_residues, (protein within 4 of ligand)
show sticks, interacting_residues

# 2. 计算相互作用
distance hydrogen_bonds, (ligand and (name O+N)), (interacting_residues and (name O+N)), cutoff=3.2

# 3. 视觉优化
set ray_shadows, 0
label interacting_residues, "%s%s" % (resn, resi)
set label_font, 12
set label_size, 20

# 4. 渲染和保存
ray 1000, 1000
png "closeup_view.png", dpi=300
```
# 5. 其他实用操作
## 5.1 文件操作
```bash
# 保存选中的结构
save "output.pdb", selection_name

# 保存会话
save "session.pse"
```
## 5.2 高级显示设置
```bash
# X-ray效果设置
Plugin - Lighting Settings - Xray

# Draw/Ray设置
set ray_trace_mode, 1
set ray_shadows, 0
```