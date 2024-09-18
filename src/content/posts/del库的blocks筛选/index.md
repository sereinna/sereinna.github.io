---
title: 根据blocks对del库筛选
published: 2024-09-18
description: "一种针对由多种blocks组成的del库的筛选，从中提取一个小的库"
image: "./121643947_p0.jpg"
tags: [Linux,del,python]
category: "代码"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 del数据库 `del`
:::

## 具体代码
```python 
import pandas as pd
import itertools
import random
from tqdm import tqdm

def select_rows(data_path, num_samples, output_to_csv=True):
    """从CSV中根据CodeA和CodeB的随机选择的值，以及CodeC为这些值加197来筛选行。"""
    data = pd.read_csv(data_path) 

    # 创建排除45和196的范围  
    excluded_values = {45, 196}  
    available_values = [x for x in range(1, 197) if x not in excluded_values]  

    # 从可用值中随机选择样本  
    values = random.sample(available_values, num_samples) 
    
    # 添加197到colC_values列表中
    colC_values = values + [197]
    print(colC_values)
    
    combinations = list(itertools.product(values, values, colC_values))
    selected_rows = []
    with tqdm(total=len(combinations)) as pbar:
        for combination in combinations:
            col1_value, col2_value, col3_value = combination
            selected_row = data[(data['Code_A'] == col1_value) &
                                (data['Code_B'] == col2_value) &
                                (data['Code_C'] == col3_value)]
            if not selected_row.empty:
                selected_rows.append(selected_row)
            pbar.update(1)
    
    if selected_rows:
        all_datasets = pd.concat(selected_rows, ignore_index=True)
    
        # 如果存在，重命名列
        if 'Exp-B01' in all_datasets.columns:
            all_datasets.rename(columns={'Exp-B01': 'Pre'}, inplace=True)
        #if 'OA' in all_datasets.columns:
            #all_datasets.rename(columns={'OA': 'Post'}, inplace=True)
        
        # 输出到CSV文件
        if output_to_csv:
            output_filename = '清洗del数据集_127500.csv'
            all_datasets.to_csv(output_filename, index=False)
            print(f"Data saved to {output_filename}")
    else:
        print("No rows found matching the criteria.")

# 使用50个样本调用函数
select_rows('清洗del数据集.csv', 50, output_to_csv=True)

```