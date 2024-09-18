---
title: g16安装脚本及集群部署
published: 2024-09-18
description: "如何在linux系统安装Gaussian高斯g16并在超算等多集群进行启动部署计算"
image: "./117158168_p0.jpg"
tags: [Gaussian,linux]
category: "代码"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 高斯安装 `Gaussian(g16)`
相关文章:[高斯安装](/posts/g16高斯安装/)[高斯运行代码](/posts/g16高斯运行代码/)
:::

## 运行安装脚本
准备文件
![g16安装准备文件](./1726664764890.jpg)
运行代码
```bash
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

#source /public/software/compiler/intel/intel-compiler-2017.5.239/bin64/compilervars.sh intel64
#source /public/software/compiler/intel/intel-compiler-2017.5.239/mkl/bin/mklvars.sh
#source /public/software/profile.d/mpi_openmpi-intel-2.1.2.sh
#source /public/software/mpi/intelmpi/2017.4.239/bin64/mpivars.sh intel64


# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/public/home/LZUZhoupp/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/public/home/LZUZhoupp/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/public/home/LZUZhoupp/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/public/home/LZUZhoupp/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# Add Gaussian 16 environment variables
export g16root=/public/home/LZUZhoupp
export GAUSS_EXEDIR=$g16root/g16
export GAUSS_SCRDIR=$g16root/g16/scr
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$g16root/g16

# Ensure Gaussian 16 binaries are in the PATH
export PATH=$PATH:$g16root/g16

# Source the Gaussian 16 profile if it exists
if [ -f "$g16root/g16/bsd/g16.profile" ]; then
    source $g16root/g16/bsd/g16.profile
fi
```

## 进行部署计算


```bash
#!/bin/bash

# 检查参数数量
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <Gaussian_input_file> <Gaussian_output_file> <Log_file>"
    exit 1
fi

# 获取输入文件名、输出文件名和日志文件名
input_file="$1"
output_file="$2"
log_file="$3"

# 检测系统上可用的 CPU 核心数
num_cores=$(grep -c ^processor /proc/cpuinfo)

# 输出一些运行信息到日志文件
echo "Running Gaussian with $num_cores cores." > "$log_file"
echo "Input file: $input_file" >> "$log_file"
echo "Output file: $output_file" >> "$log_file"
echo "Log file: $log_file" >> "$log_file"

# 创建一个临时的输入文件，包含正确的 %NProcShared 参数
temp_input_file="temp_${input_file}"
echo "%NProcShared=${num_cores}" > "$temp_input_file"
cat "$input_file" >> "$temp_input_file"

# 运行 Gaussian 并将输出重定向到指定的输出文件
g16 < "$temp_input_file" > "$output_file" 2>> "$log_file"

# 记录完成信息到日志文件
echo "Gaussian run completed." >> "$log_file"

# 清理临时文件
#rm "$temp_input_file"
```
