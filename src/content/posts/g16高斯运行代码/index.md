---
title: g16运行脚本
published: 2024-09-18
description: "如何使用Gaussian高斯g16去计算一个分子的xx"
image: "./121711768_p0.png"
tags: [Gaussian,python]
category: "代码"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 高斯安装 `Gaussian(g16)`
相关文章:[高斯安装](/posts/g16高斯安装/) | [高斯运行代码](/posts/g16高斯运行代码/)
:::

## 运行脚本
```python 
import subprocess
import os
import shutil

files = ["hth-6-td-TICT-OPT.gjf"]

for file_path in files:

    base_dir = os.path.dirname(file_path)
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_dir = os.path.join(base_dir, f"{base_name}_results")
    os.makedirs(output_dir, exist_ok=True)
    shutil.copy(file_path, output_dir)
    file_path_in_new_dir = os.path.join(output_dir, os.path.basename(file_path))

    # 设置环境变量最大内存
    os.environ['GAUSS_PDEF'] = '48'
    os.environ['GAUSS_MEMDEF'] = '28GB'

    # 构建Gaussian命令
    out_file_path = os.path.join(output_dir, f"{base_name}.out")
    g16_command = f"g16 < {file_path_in_new_dir} > {out_file_path}"
    print(f"Starting Gaussian calculation for {base_name}...")
    subprocess.run(g16_command, shell=True, check=True)


    # 构建formchk命令
    #chk_file = os.path.join(output_dir, f"{base_name}.chk")
    #fchk_file = os.path.join(output_dir, f"{base_name}.fchk")
    #formchk_command = f"formchk {chk_file} {fchk_file}"
    #print(f"Converting .chk to .fchk for {base_name}...")

    print(f"Completed processing for {base_name}.")

print("All files processed.")

```