---
title: 分子对接
published: 2024-09-18
description: "如何使用smina批量对蛋白质-配体分子进行对接"
image: "./122503112_p0.png"
tags: [python,docking,smina]
category: "代码"
draft: false
lang: ""
---

## 前言

:::tip
准备工作  配体文件`sdf`，蛋白文件`pdb`\

:::



## 对接代码
```python 
import subprocess
import os
from tqdm import tqdm
from multiprocessing import Pool
from openbabel import pybel


def run_smina_docking(params):
    input_path, output_sdf, protein, ligand, box_size, num_results, num_cores = params
    # Check if the output file already exists and check the number of conformations
    if os.path.exists(output_sdf) and os.path.getsize(output_sdf) > 0:
        try:
            # Load the SDF file and count the conformations
            conformations = [mol for mol in pybel.readfile('sdf', output_sdf)]
            if len(conformations) >= num_results:
                return
        except Exception as e:
            print(f"Error reading {output_sdf}: {e}")
            # You might want to handle this differently based on your needs

    command_sdf = [
        'smina', '--receptor', protein, '--ligand', input_path,
        '--autobox_ligand', ligand, '--autobox_add', str(box_size),
        '--num_modes', str(num_results), '--cpu', str(num_cores),
        '--out', output_sdf
    ]
    subprocess.run(command_sdf, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
def run_docking(input_folder, output_folder, protein_name, box_size, num_results, num_cores, batch_size=2000, start_idx=None, end_idx=None):
    input_folder = f'{input_folder}'
    output_folder = f'{output_folder}_{protein_name}'
    protein = f'pdb/{protein_name}_protein.pdb'
    ligand = f'pdb/{protein_name}_ligand.sdf'

    os.makedirs(output_folder, exist_ok=True)

    sdf_files = sorted([f for f in os.listdir(input_folder) if f.endswith('.sdf')], key=lambda x: int(os.path.splitext(x)[0]))

    if start_idx is not None and end_idx is not None:
        sdf_files = sdf_files[start_idx-1:end_idx]
    elif start_idx is not None:
        sdf_files = sdf_files[start_idx-1:]
    elif end_idx is not None:
        sdf_files = sdf_files[:end_idx]

    total_files = len(sdf_files)

    for batch_start in range(0, total_files, batch_size):
        batch_end = min(batch_start + batch_size, total_files)
        batch_files = sdf_files[batch_start:batch_end]

        params_list = []
        for sdf_file in batch_files:
            input_path = os.path.join(input_folder, sdf_file)
            output_base = os.path.splitext(sdf_file)[0]
            output_sdf = os.path.join(output_folder, f"{output_base}_out.sdf")
            params_list.append((input_path, output_sdf, protein, ligand, box_size, num_results, num_cores))

        with Pool(processes=num_cores) as pool:
            for _ in tqdm(pool.imap_unordered(run_smina_docking, params_list), total=len(batch_files), desc="Docking progress", unit="files"):
                pass

        print(f"Batch {batch_start // batch_size + 1} completed. Processed files {batch_start+1} to {batch_end}.")

    print("Docking process completed.")

# Example usage
protein_name = '3p3h'
input_folder = '清洗数据集_0720_input'
output_folder = 'dock/output_clean'
box_size = 25
num_results = 20
num_cores = 40
batch_size = 5000  # Set the batch size to 2000

# Run docking for ligands in the range 1 to 50000
run_docking(input_folder, output_folder, protein_name, box_size, num_results, num_cores, batch_size=batch_size, start_idx=1, end_idx=127500)
```