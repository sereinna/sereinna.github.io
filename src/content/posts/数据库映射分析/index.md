---
title: 基于Tanimoto核的t-SNE算法
published: 2024-09-23
description: "使用基于Tanimoto核的t-SNE算法对蛋白分子数据库可视化分子空间的高维特征"
image: "./120005527_p0_master1200.jpg"
tags: [数据集,分子,可视化,python]
category: "代码"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 多个蛋白/分子数据集 `DrugBank or BIOSNAP or BindingDB`，
[相关文献](https://mp.weixin.qq.com/s/mZbTyZ-c6v7TQMIAZ3Y4cA)


## 基于Tanimoto核的t-SNE算法

基于分子特征使用t-SNE算法对分子空间进行二维表示是一种有效的可视化方法。化学空间被认为极其庞大，训练数据集中对该空间的高效采样预计将对预测模型的泛化性能产生重大影响。为可视化分子空间的高维特征，研究者常使用基于Tanimoto核的t-SNE算法，将其映射到2D空间中进行分析。

### t-SNE算法

t-SNE（t-distributed Stochastic Neighbor Embedding）是一种非线性降维算法，特别适合用于高维数据的可视化。它能够保持数据点之间的局部关系，同时揭示数据的全局结构。

### 分子空间

化学空间是指所有可能的化学结构的集合。这个空间是高维的，因为每个分子可以用多个特征来描述（如原子类型、化学键、官能团等）。

### Tanimoto核

Tanimoto系数（也称为Jaccard系数）是衡量两个集合相似度的指标。在化学信息学中，它常用于比较分子指纹的相似性。Tanimoto核是基于这个系数的核函数，用于在高维空间中计算分子间的相似度。

### 高效采样的重要性

由于化学空间极其庞大，无法穷尽所有可能的分子。因此，训练数据集中的分子应该能够代表整个化学空间的多样性，这对于模型的泛化性能至关重要。

### Python代码实现

以下是一个基于Python的代码示例，展示如何使用t-SNE算法对分子特征进行降维和可视化：

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

# 1. 准备分子数据
smiles_list = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",  # Testosterone
    "CN1CCC[C@H]1C2=CN=CC=C2",  # Nicotine
    "CC(C)(C)NCC(O)C1=CC(=C(C=C1)O)CO"  # Salbutamol
]

# 2. 计算分子指纹
def calculate_fingerprints(smiles_list):
    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]
    return np.array(fps)

fingerprints = calculate_fingerprints(smiles_list)

# 3. 定义Tanimoto核函数
def tanimoto_kernel(X):
    return 1 - pdist(X, metric="jaccard")

# 4. 使用t-SNE进行降维
tsne = TSNE(n_components=2, metric="precomputed", random_state=42)
tanimoto_dist = squareform(tanimoto_kernel(fingerprints))
tsne_results = tsne.fit_transform(1 - tanimoto_dist)

# 5. 可视化结果
plt.figure(figsize=(10, 8))
plt.scatter(tsne_results[:, 0], tsne_results[:, 1])
for i, smiles in enumerate(smiles_list):
    plt.annotate(Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=False), 
                 (tsne_results[i, 0], tsne_results[i, 1]))
plt.title("t-SNE visualization of molecular space")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.show()
```