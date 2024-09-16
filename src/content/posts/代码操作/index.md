---
title: 操作代码记录
published: 2024-09-16
description: "Linux系统，Del数据，docking"
image: "./1.jpg"
tags: [Linux,del,docking,代码]
category: "代码"
draft: false
lang: ""
---

## Linux
```bash
#多显卡运行
CUDA_VISIBLE_DEVICES = "0,1,2,3" 
```

```bash
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
del信息
P00918 · CAH2_HUMAN      pdb_id: 3p3h 5doh
Q16790 · CAH9_HUMAN     pdb_id: 2hkf 5fl4
O43570 · CAH12_HUMAN   pdb_id: 4kp5 4ht2


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


## 高斯运行
```bash
#高斯运行
qsub g16.pbs
运行
```

## 高斯计算
```markdown  
```python 
import subprocess
import os
import shutil

files = ["hth7-td-san-1.gjf"]

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


## win激活.reg
```markdown  
```python 
Windows Registry Editor Version 5.00

[HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\.NETFramework\v2.0.50727]
"SystemDefaultTlsVersions"=dword:00000001
"SchUseStrongCrypto"=dword:00000001

[HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\.NETFramework\v4.0.30319]
"SystemDefaultTlsVersions"=dword:00000001
"SchUseStrongCrypto"=dword:00000001

[HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\.NETFramework\v2.0.50727]
"SystemDefaultTlsVersions"=dword:00000001
"SchUseStrongCrypto"=dword:00000001

[HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\.NETFramework\v4.0.30319]
"SystemDefaultTlsVersions"=dword:00000001
"SchUseStrongCrypto"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols]

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\SSL 2.0]

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\SSL 2.0\Client]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\SSL 2.0\Server]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\SSL 3.0]

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\SSL 3.0\Client]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\SSL 3.0\Server]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.0]

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.0\Client]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.0\Server]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.1]

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.1\Client]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.1\Server]
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.2]
@=""

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.2\Client]
@=""
"DisabledByDefault"=dword:00000000
"Enabled"=dword:00000001

[HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\SecurityProviders\SCHANNEL\Protocols\TLS 1.2\Server]
@=""
"Enabled"=dword:00000001
"DisabledByDefault"=dword:00000000

```

## Chembl数据筛选
```markdown  
```python 
import csv
from tqdm import tqdm
from chembl_webresource_client.new_client import new_client

# 连接到 ChEMBL 数据库
target = new_client.target
activity = new_client.activity
molecule = new_client.molecule

# 使用 UniProt ID 搜索目标蛋白
target_query = target.search('O43570')
targets = [t for t in target_query if t['target_components'][0]['accession'] == 'O43570']
if targets:
    target_id = targets[0]['target_chembl_id']

    # 获取与目标蛋白相关的活性数据
    activities = activity.filter(target_chembl_id=target_id).filter(standard_type__in=['Ki', 'IC50'])

    # 筛选报告了 Ki 和 IC50 值的化合物
    compounds = [a for a in activities if a['standard_value'] is not None]

    # 准备保存到 CSV 文件的数据
    compounds_data = []

    # 使用 tqdm 创建进度条
    for compound in tqdm(compounds, desc="处理化合物"):
        # 获取化合物的 SMILES 字符串
        compound_record = molecule.get(compound['molecule_chembl_id'])
        smiles = compound_record['molecule_structures']['canonical_smiles'] if compound_record['molecule_structures'] else 'N/A'

        # 添加到数据列表
        compounds_data.append([
            compound['molecule_chembl_id'],
            compound['standard_type'],
            compound['standard_value'],
            compound['standard_units'],
            smiles
        ])

    # 将结果保存到 CSV 文件
    with open('compounds_activities.csv', 'w', newline='', encoding='utf-8-sig') as file:
        writer = csv.writer(file)
        # 写入表头
        writer.writerow(['化合物 ChEMBL ID', '活性类型', '活性值', '单位', 'SMILES'])

        # 写入数据
        writer.writerows(compounds_data)

    print("数据已保存到 compounds_activities.csv")
    ```