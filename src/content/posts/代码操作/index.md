---
title: 操作代码记录
published: 2024-09-16
description: "Linux系统，Del数据，docking"
image: ""
tags: [Linux,del,docking]
category: "我如何使用gpt获得一堆代码并记录"
draft: false
lang: ""
---

## Linux
```bash
#多显卡运行
CUDA_VISIBLE_DEVICES = "0,1,2,3" 

#清华源配置
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch-lts: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  deepmodeling: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/

```


## conda


## Del
```bash
#del信息
P00918 · CAH2_HUMAN      pdb_id: 3p3h 5doh
Q16790 · CAH9_HUMAN     pdb_id: 2hkf 5fl4
O43570 · CAH12_HUMAN   pdb_id: 4kp5 4ht2

#

```


## Docking
```bash
#dock环境安装
sudo apt update
sudo apt install openssh-server
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda create -n dock

conda activate dock
conda install -c conda-forge smina
conda install tqdm
#

```


## 高斯
```bash
#高斯运行
qsub g16.pbs
运行
```