---
title: 这里写标题
published: 2024-09-17
description: "这是一个撰写模板"
image: "./1.jpg"
tags: [Linux,del,这里写标签]
category: "这里写分类"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 Cloudflare `域名`，NPM `反向代理`，SSL `证书`\
相关文章:[反向代理](/posts/npm-install/)|[SSL 证书](/posts/acme/)

:::

## 
```python 
#c:\Users\86173\Documents\WeChat Files\wxid_azbxtur2sd4x22\FileStorage\Temp\5490076f7f985a48d73514f181bbd4e.png
#箱线图
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# 初始数据：方法和指标
methods = ['Dos-del', 'Deldenoiser', 'DEL-MLE', 'DEL-MAP', 'MPDF']
metrics = ['AUROC', 'AUPRC', 'Recall', 'F1', 'Prec']

# 每个方法的均值和误差
means = [
    [0.923, 0.727, 0.877, 0.832, 0.830],  # Dos-del
    [0.840, 0.150, 0.750, 0.269, 0.164],  # Deldenoiser
    [0.871, 0.766, 0.741, 0.849, 1.000],  # DEL-MLE
    [1.000, 1.000, 1.000, 1.000, 1.000],  # DEL-MAP
    [1.000, 1.000, 1.000, 1.000, 1.000],  # MPDF
]

errors = [
    [0.113, 0.193, 0.215, 0.170, 0.141],  # Dos-del
    [0.070, 0.074, 0.144, 0.083, 0.081],  # Deldenoiser
    [0.051, 0.088, 0.102, 0.065, 0.000],  # DEL-MLE
    [0.000, 0.000, 0.000, 0.000, 0.000],  # DEL-MAP
    [0.000, 0.000, 0.000, 0.000, 0.000],  # MPDF
]

# 设置每个结果的颜色
colors = ['#845ec2', '#4b4453', '#b0a8b9', '#00896f', '#00c0a3']

# 创建箱线图
fig, ax = plt.subplots(figsize=(10, 6))

# 计算每种方法在x轴上的位置
num_metrics = len(metrics)
positions = np.arange(len(methods)) * (num_metrics + 1)  # 为每种方法分配x轴上的空间

# 为每种方法生成箱线图，每个方法内部的五个结果使用不同颜色
for i, method_means in enumerate(means):
    boxplot_data = []
    for j in range(num_metrics):
        lower_bound = method_means[j] - errors[i][j]
        upper_bound = method_means[j] + errors[i][j]
        # 创建数据
        boxplot_data.append([lower_bound, method_means[j], upper_bound])
    
    # 绘制箱线图，每个指标分配不同的颜色
    bp = ax.boxplot(boxplot_data, positions=positions[i] + np.arange(num_metrics), patch_artist=True, widths=0.6)
    
    # 为每个箱体设置颜色
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    # 修改中值线的颜色和长度
    for median in bp['medians']:
        median.set_color('#334a52')  # 修改中值线的颜色
        median.set_linewidth(1.5)  # 保持中值线的粗细
        median.set_xdata([median.get_xdata()[0] + 0.05, median.get_xdata()[1] - 0])  # 确保中值线不会突出

# 设置x轴标签
ax.set_xticks([p + (num_metrics / 2) - 0.5 for p in positions])  # 调整x轴标签的位置
ax.set_xticklabels(methods, fontsize=12, fontweight='bold')

# 设置y轴标签
ax.set_ylabel('Scores', fontsize=12, fontweight='bold')

# 设置标题
ax.set_title('BB10', fontsize=14, fontweight='bold')

# 设置轴刻度标签加粗
ax.tick_params(axis='both', labelsize=12, width=2)

# 保存箱线图为SVG文件
fig.tight_layout()
plt.savefig('/home/tim/hmt/del_picture/误差图/database/箱线图.svg', format='svg')

# 显示箱线图
plt.show()

# 创建图例图
fig_legend, ax_legend = plt.subplots(figsize=(6, 3))

# 添加填充的长方形作为图例
handles = []
for color, metric in zip(colors, metrics):
    rect = Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black')  # 创建有边框的长方形
    handles.append(rect)

# 显示图例
ax_legend.legend(handles=handles, labels=metrics, title="Metrics", title_fontsize=12, fontsize=10, loc="center", frameon=True, edgecolor='black')
ax_legend.axis('off')  # 隐藏坐标轴

# 保存图例为SVG文件
fig_legend.tight_layout()
plt.savefig('/home/tim/hmt/del_picture/误差图/database/箱线图图例.svg', format='svg')

# 显示图例
plt.show()
```


图片参考 `T-SNE`
![部署图片](./1.jpg)
