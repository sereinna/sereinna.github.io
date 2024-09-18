---
title: sdf文件生成
published: 2024-09-17
description: "如何批量根据smiles生成sdf 3d文件"
image: "./8.jpg"
tags: [python,docking]
category: "代码"
draft: false
lang: ""
---

## 前言

:::tip
准备工作  配体`smiles`，蛋白 `smiles`\

:::

## 代码
```python 
import pandas as pd
from openbabel import pybel
import os
from tqdm import tqdm

# 读取csv文件
data = pd.read_csv('filtered_Q16790.csv')

# 获取总行数
total_rows = len(data)

# 创建目标目录，如果不存在的话
output_dir = '/home/zjlab/dock/8w_yanzhengji_dock/Q16790/docking_input'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 遍历每一行，显示进度条
for index, row in tqdm(data.iterrows(), total=total_rows):
    smiles = row['smiles']
    index_num = row['index']  # 假设这里的'index'是你想要的列名
    output_file = os.path.join(output_dir, f'docked_{index_num}.sdf')
    
    # 检查目标文件是否已存在
    if os.path.exists(output_file):
        continue

    # 生成sdf文件
    mol = pybel.readstring('smi', smiles)
    mol.addh()
    mol.make3D()
    mol.write('sdf', output_file)
```