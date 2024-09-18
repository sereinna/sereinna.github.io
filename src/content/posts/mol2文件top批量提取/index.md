---
title: 批量提取top_mol2
published: 2024-09-18
description: "如何使用python代码批量提取对接后的mol2文件中的top1"
image: "./121653445_p0.jpg"
tags: [python,docking]
category: "代码"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 分子文件 `mol2`
:::

## 
```python 
import os
import re

def extract_molecule(input_file_path, output_file_path):
    with open(input_file_path, 'r') as input_file:
        mol_lines = []
        count = 0
        for line in input_file:
            if line.startswith('@<TRIPOS>MOLECULE') and count == 1:
                with open(output_file_path, 'w') as output_file:
                    for mol_line in mol_lines:
                        output_file.write(mol_line)
                    break
            elif line.startswith('@<TRIPOS>MOLECULE'):
                count += 1
            mol_lines.append(line)


def process_sdf_files(root_folder):
    for folder_name in os.listdir(root_folder):
        docking_input_folder = os.path.join(root_folder, folder_name, 'ligand_smina_poses')
        if not os.path.exists(docking_input_folder):
            continue
        ligand_smina_poses_folder = os.path.join(root_folder, folder_name, 'ligand_smina_poses')
        os.makedirs(ligand_smina_poses_folder, exist_ok=True)
        for file_name in os.listdir(docking_input_folder):
            if re.match(r'^\d+\.mol2$', file_name):
                input_file_path = os.path.join(docking_input_folder, file_name)
                output_file_path = os.path.join(ligand_smina_poses_folder, file_name[:-5] + '_1.mol2')
                extract_molecule(input_file_path, output_file_path)
                print(f'Processed {input_file_path}.')
                print(f'Processed {output_file_path}.')

root_folder = 'only_chembl_3'
process_sdf_files(root_folder)

```