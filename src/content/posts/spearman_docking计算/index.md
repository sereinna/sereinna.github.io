---
title: 计算对接结果的spearman常数
published: 2024-09-17
description: "python的三种spearman常数计算过程"
image: "./3.jpg"
tags: [spearman,python,docking]
category: "代码"
draft: false
lang: ""
---

## 如何计算一个对接结果的benzene的spearman常数
```python 
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm

# 读取 CSV 文件
# df_eval = pd.read_csv("df_eval_data_P00918.csv")
df_eval = pd.read_csv("/home/tim/hmt/zj/验证集/our/P00918/df_eval_data_Q16790_5fl4.csv")

# 过滤掉 'Ki (nM)' 列中包含 NaN 值的行
df_eval = df_eval.dropna(subset=['Ki (nM)'])

# 定义检测 benzene sulfonamide 的函数
def has_benzene_sulfonamide(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol.HasSubstructMatch(Chem.MolFromSmarts('c1ccccc1S(=O)(=O)N'))

# 计算每个分子的分子量和 benzene sulfonamide 存在与否
mws = []
benzene_sulfonamide_presence = []
for smi in tqdm(df_eval["smiles"], desc="计算分子量和检测 benzene sulfonamide"):
    mol = Chem.MolFromSmiles(smi)
    mws.append(Descriptors.MolWt(mol))
    benzene_sulfonamide_presence.append(int(has_benzene_sulfonamide(smi)))

mws = np.array(mws)
benzene_sulfonamide_presence = np.array(benzene_sulfonamide_presence)

# 获取 Ki 值
ki_values = df_eval['Ki (nM)'].values

# 计算全数据集的 Spearman 相关系数
spearman_full, _ = spearmanr(benzene_sulfonamide_presence, ki_values)
print("Spearman Ki (full) ↓:", spearman_full)

# 计算分子量的第 10 百分位数和第 90 百分位数
lower_bound = np.percentile(mws, 10)
upper_bound = np.percentile(mws, 90)

# 根据计算出的分子量范围筛选分子
valid_idxs = np.logical_and(mws > lower_bound, mws < upper_bound)

# 筛选符合条件的 benzene sulfonamide 存在与否和对应的 Ki 值
filtered_benzene_sulfonamide = benzene_sulfonamide_presence[valid_idxs]
filtered_ki_values = ki_values[valid_idxs]

# 过滤掉筛选后的 Ki 值中的 NaN
filtered_valid_idxs = ~np.isnan(filtered_ki_values)
filtered_benzene_sulfonamide = filtered_benzene_sulfonamide[filtered_valid_idxs]
filtered_ki_values = filtered_ki_values[filtered_valid_idxs]

# 计算筛选后数据集的 Spearman 相关系数
if len(filtered_benzene_sulfonamide) > 1 and len(np.unique(filtered_benzene_sulfonamide)) > 1 and len(np.unique(filtered_ki_values)) > 1:
    spearman_filtered, _ = spearmanr(filtered_benzene_sulfonamide, filtered_ki_values)
    print(f"Spearman Ki ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓:", spearman_filtered)
else:
    print(f"Spearman Ki ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓: nan (due to insufficient variation or sample size)")

#print("Unique values in filtered Ki values:", np.unique(filtered_ki_values))

```


## 如何计算一个对接结果的weight的spearman常数
```python 
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm  # 引入 tqdm

# 读取 CSV 文件，假设 CSV 中有 'smiles' 和 'Ki (nM)' 两列
df_eval = pd.read_csv("/home/tim/hmt/zj/验证集/our/P00918/df_eval_data_Q16790_5fl4.csv")

# 过滤掉 'Ki (nM)' 和 'smiles' 列中包含 NaN 的行
df_eval = df_eval.dropna(subset=['Ki (nM)', 'smiles'])

# 计算每个分子的分子量
mws = np.array([Descriptors.MolWt(Chem.MolFromSmiles(smi)) for smi in tqdm(df_eval["smiles"], desc="计算分子量")])

# 获取 Ki 值
ki_values = df_eval['Ki (nM)'].values

# 计算 Spearman Ki (full) ↓，即针对所有分子的 Spearman 相关系数
spearman_full, _ = spearmanr(mws, ki_values)
print("Spearman Ki (full) ↓:", spearman_full)

# 计算分子量的第 10 百分位数和第 90 百分位数
lower_bound = np.percentile(mws, 10)
upper_bound = np.percentile(mws, 90)

# 根据计算出的分子量范围筛选分子
valid_idxs = np.logical_and(mws > lower_bound, mws < upper_bound)

# 筛选符合条件的分子量和对应的 Ki 值
filtered_mws = mws[valid_idxs]
filtered_ki_values = ki_values[valid_idxs]

# 过滤掉筛选后的 Ki 值中的 NaN
filtered_valid_idxs = ~np.isnan(filtered_ki_values)
filtered_mws = filtered_mws[filtered_valid_idxs]
filtered_ki_values = filtered_ki_values[filtered_valid_idxs]

# 计算 Spearman Ki (10th < MW < 90th) ↓，即针对分子量在第 10 和第 90 百分位之间的 Spearman 相关系数
if len(filtered_mws) > 1 and len(np.unique(filtered_mws)) > 1 and len(np.unique(filtered_ki_values)) > 1:
    spearman_filtered, _ = spearmanr(filtered_mws, filtered_ki_values)
    print(f"Spearman Ki ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓:", spearman_filtered)
else:
    print(f"Spearman Ki ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓: nan (due to insufficient variation or sample size)")

```


## 如何计算一个对接结果的weight的spearman常数
```python 
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm  # 引入 tqdm

# 读取 CSV 文件，假设 CSV 中有 'smiles' 和 'Ki (nM)' 两列
df_eval = pd.read_csv("/home/tim/hmt/zj/验证集/our/P00918/df_eval_data_Q16790_5fl4.csv")

# 过滤掉 'Ki (nM)' 和 'smiles' 列中包含 NaN 的行
df_eval = df_eval.dropna(subset=['Ki (nM)', 'smiles'])

# 计算每个分子的分子量
mws = np.array([Descriptors.MolWt(Chem.MolFromSmiles(smi)) for smi in tqdm(df_eval["smiles"], desc="计算分子量")])

# 获取 Ki 值
ki_values = df_eval['Ki (nM)'].values

# 计算 Spearman Ki (full) ↓，即针对所有分子的 Spearman 相关系数
spearman_full, _ = spearmanr(mws, ki_values)
print("Spearman Ki (full) ↓:", spearman_full)

# 计算分子量的第 10 百分位数和第 90 百分位数
lower_bound = np.percentile(mws, 10)
upper_bound = np.percentile(mws, 90)

# 根据计算出的分子量范围筛选分子
valid_idxs = np.logical_and(mws > lower_bound, mws < upper_bound)

# 筛选符合条件的分子量和对应的 Ki 值
filtered_mws = mws[valid_idxs]
filtered_ki_values = ki_values[valid_idxs]

# 过滤掉筛选后的 Ki 值中的 NaN
filtered_valid_idxs = ~np.isnan(filtered_ki_values)
filtered_mws = filtered_mws[filtered_valid_idxs]
filtered_ki_values = filtered_ki_values[filtered_valid_idxs]

# 计算 Spearman Ki (10th < MW < 90th) ↓，即针对分子量在第 10 和第 90 百分位之间的 Spearman 相关系数
if len(filtered_mws) > 1 and len(np.unique(filtered_mws)) > 1 and len(np.unique(filtered_ki_values)) > 1:
    spearman_filtered, _ = spearmanr(filtered_mws, filtered_ki_values)
    print(f"Spearman Ki ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓:", spearman_filtered)
else:
    print(f"Spearman Ki ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓: nan (due to insufficient variation or sample size)")

```



## 如何统计一个对接结果的value
```python 
import pandas as pd
import os

# 路径定义
csv_file_path = '/home/tim/hmt/zj/验证集/our/P00918/df_eval_data_Q16790_5fl4.csv'
#csv_file_path = '/home/tim/hmt/zj/验证集/our/P00918/df_eval_data_Q16790_5fl6.csv'

log_dir_path = '/home/tim/hmt/zj/验证集/our/Q16790/5fl4'
#log_dir_path = '/home/tim/hmt/zj/验证集/our/P00918/5doh'

# 读取CSV文件
df = pd.read_csv(csv_file_path)

# 为每个affinity值添加新列
for i in range(1, 10):
    df[f'top_docking_affinity_{i}'] = None

# 遍历DataFrame的行
for index, row in df.iterrows():
    num = row['num']
    log_file_name = f"docked_{num}.log"
    log_file_path = os.path.join(log_dir_path, log_file_name)
    
    # 检查文件是否存在
    if os.path.isfile(log_file_path):
        with open(log_file_path, 'r') as file:
            pose_section = False
            mode_count = 1
            for line in file:
                if "mode |   affinity | dist from best mode" in line:
                    pose_section = True  # 开始pose部分
                    continue
                if pose_section and line.strip():
                    parts = line.split()
                    if parts[0].isdigit() and mode_count <= 9:
                        try:
                            affinity = float(parts[1])  # 提取affinity值
                            df.at[index, f'top_docking_affinity_{mode_count}'] = affinity
                            mode_count += 1
                        except ValueError:
                            continue

# 保存更新后的CSV文件
df.to_csv(csv_file_path, index=False)

```

## 如何计算一个对接结果的value的spearman常数
```python 
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm

# 读取CSV文件
df_eval = pd.read_csv("/home/tim/hmt/zj/验证集/our/P00918/df_eval_data_Q16790_5fl4.csv")

# 过滤掉包含NaN的行
df_eval = df_eval.dropna(subset=['Ki (nM)', 'smiles', 'top_docking_affinity_1'])

# 将top_docking_affinity转换为正数
df_eval['top_docking_affinity_1'] = df_eval['top_docking_affinity_1'].abs()

# 添加每个affinity值的列
affinity_columns = [f'top_docking_affinity_{i}' for i in range(1, 10)]
df_eval[affinity_columns] = df_eval[affinity_columns].abs()  # 确保所有亲和力值为正

# 计算平均值
df_eval['average_affinity'] = df_eval[affinity_columns].mean(axis=1)

# 计算最大值和最小值的中位数
df_eval['median_of_extremes'] = df_eval[affinity_columns].apply(lambda x: np.median([x.max(), x.min()]), axis=1)

# 计算每个分子的分子量
df_eval['molecular_weight'] = [Descriptors.MolWt(Chem.MolFromSmiles(smi)) for smi in tqdm(df_eval["smiles"], desc="计算分子量")]

# 获取Ki值和其他affinity列的值
ki_values = df_eval['Ki (nM)'].values
top_docking_affinities = df_eval['top_docking_affinity_1'].values
average_affinities = df_eval['average_affinity'].values
median_extremes = df_eval['median_of_extremes'].values

# 计算Spearman相关系数
spearman_original, _ = spearmanr(top_docking_affinities, ki_values)
spearman_average, _ = spearmanr(average_affinities, ki_values)
spearman_median_extremes, _ = spearmanr(median_extremes, ki_values)

# 计算分子量的第10和90百分位
lower_bound = np.percentile(df_eval['molecular_weight'], 10)
upper_bound = np.percentile(df_eval['molecular_weight'], 90)

# 筛选符合条件的分子量范围内的数据
filtered_df = df_eval[(df_eval['molecular_weight'] > lower_bound) & (df_eval['molecular_weight'] < upper_bound)]

# 计算筛选后数据的Spearman相关系数
filtered_spearman_original, _ = spearmanr(filtered_df['top_docking_affinity_1'], filtered_df['Ki (nM)'])
filtered_spearman_average, _ = spearmanr(filtered_df['average_affinity'], filtered_df['Ki (nM)'])
filtered_spearman_median_extremes, _ = spearmanr(filtered_df['median_of_extremes'], filtered_df['Ki (nM)'])

# 打印所有Spearman相关系数结果
print("Spearman Ki (full) ↓:", spearman_original)
print("Spearman Ki (average) ↓:", spearman_average)
print("Spearman Ki (median of extremes) ↓:", spearman_median_extremes)
print(f"Spearman Ki filtered ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓ - Original:", filtered_spearman_original)
print(f"Spearman Ki filtered ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓ - Average:", filtered_spearman_average)
print(f"Spearman Ki filtered ({lower_bound:.2f} < MW < {upper_bound:.2f}) ↓ - Median of Extremes:", filtered_spearman_median_extremes)
```