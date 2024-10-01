---
title: smiles分子图批量绘制
published: 2024-10-01
description: "如何根据几个选定的颜色绘制还算好看的smlies图，并突出苯磺酰胺基团"
image: "./122156399_p0.jpg"
tags: [python,分子,图]
category: "画图"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 分子 `smiles`

:::

## 
```python 
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import os
from rdkit.Chem import inchi

# 原子高亮和颜色设置（可以在此快速设置）
#HIGHLIGHT_ATOMS = ['N', 'S', 'O', 'F', 'Cl', 'Br', 'I']  
HIGHLIGHT_ATOMS = [] 
ATOM_COLORS = {
    'N': (1.0, 0.7569, 0.6235),  # 氮原子的颜色 FFC19F
    'S': (0.9725, 0.7059, 0.8627),  # 硫原子的颜色 F8B4DC
    'O': (0.0000, 0.7882, 0.7176),   # 氧原子的颜色 #008578 深青绿色
    'F': (0.9725, 0.7059, 0.8627),   # 氟原子的颜色 F8B4DC
    'Cl': (0.1216, 0.6824, 0.2824),  # 氯原子的颜色 #1FAD48 绿色
    'Br': (0.6510, 0.1647, 0.1647),  # 溴原子的颜色 #A62929 深红色
    'I': (0.4941, 0.0980, 0.6157),   # 碘原子的颜色 #7E199D 紫色
}

# 特殊基团颜色（例如苯磺酰胺）
BENZENESULFONAMIDE_COLOR = (0.9922, 0.5765, 0.7686)  # #fd93c4 亮粉红色
# 高亮键的颜色
HIGHLIGHT_BOND_COLOR = (1.0000, 0.7941, 0.9725)  # #ffe4f8 非常浅的粉红色

# 自动生成颜色映射
def assign_colors_by_atom_type(mol, highlight_atoms, atom_colors):
    color_map = {}
    for i, atom in enumerate(mol.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        if atom_symbol in highlight_atoms and atom_symbol in atom_colors:
            color_map[i] = atom_colors[atom_symbol]
    return color_map

# 识别苯磺酰胺基团的函数
def identify_benzenesulfonamide(mol):
    pattern = Chem.MolFromSmarts('c1ccc(S(=O)(N)=O)cc1')
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return list(matches[0])
    return []

# 辅助函数：生成唯一的文件名
def generate_unique_filename(directory, base_filename, extension):
    """
    如果文件名已存在，则在文件名后添加下标以确保唯一性。
    例如，如果 base_filename.svg 已存在，则生成 base_filename_1.svg, base_filename_2.svg 等。
    """
    filename = f"{base_filename}.{extension}"
    filepath = os.path.join(directory, filename)
    counter = 1
    while os.path.exists(filepath):
        filename = f"{base_filename}_{counter}.{extension}"
        filepath = os.path.join(directory, filename)
        counter += 1
    return filename

def draw_and_save_smiles(smiles, svg_dir, jpg_dir, idx, highlight_atoms, atom_colors):
    # 创建分子对象
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        print(f"Invalid SMILES at index {idx}: {smiles}")
        return None
    
    # 生成2D坐标
    AllChem.Compute2DCoords(mol)
    
    # 绘制设置
    draw_options = rdMolDraw2D.MolDrawOptions()
    draw_options.useBWAtomPalette()  # 关闭默认颜色
    draw_options.fixedBondLength = 30
    draw_options.bondLineWidth = 1  # 设置化学键线条的宽度为 1，值越小线条越细
    
    # 准备绘图和颜色映射
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.SetDrawOptions(draw_options)
    color_map = assign_colors_by_atom_type(mol, highlight_atoms, atom_colors)
    
    # 识别苯磺酰胺基团
    benzenesulfonamide_atoms = identify_benzenesulfonamide(mol)
    
    # 添加苯磺酰胺基团的颜色
    for atom in benzenesulfonamide_atoms:
        color_map[atom] = BENZENESULFONAMIDE_COLOR
    
    # 高亮的原子
    highlight_atoms_idx = list(set([i for i, atom in enumerate(mol.GetAtoms()) 
                                    if atom.GetSymbol() in highlight_atoms] + benzenesulfonamide_atoms))
    
    # 高亮的键
    highlight_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in benzenesulfonamide_atoms and bond.GetEndAtomIdx() in benzenesulfonamide_atoms:
            highlight_bonds.append(bond.GetIdx())
    
    # 绘制分子
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms_idx, highlightAtomColors=color_map,
                        highlightBonds=highlight_bonds, highlightBondColors={i: HIGHLIGHT_BOND_COLOR for i in highlight_bonds})
    drawer.FinishDrawing()
    svg_output = drawer.GetDrawingText()
    
    # 生成InChIKey，并提取中间部分
    inchi_key = inchi.MolToInchiKey(mol)
    base_filename = inchi_key[4:12]  # 提取中间部分
    
    # 生成唯一的SVG文件名
    unique_svg_filename = generate_unique_filename(svg_dir, base_filename, "svg")
    svg_filepath = os.path.join(svg_dir, unique_svg_filename)
    
    # 保存SVG文件
    with open(svg_filepath, "w") as svg_file:
        svg_file.write(svg_output.replace('svg:', ''))
    
    # 生成唯一的JPG文件名
    unique_jpg_filename = generate_unique_filename(jpg_dir, base_filename, "jpg")
    jpg_filepath = os.path.join(jpg_dir, unique_jpg_filename)
    
    # 保存JPG文件
    img = Draw.MolToImage(mol, size=(900, 900), highlightAtoms=highlight_atoms_idx, highlightAtomColors=color_map,
                          highlightBonds=highlight_bonds)
    img.save(jpg_filepath)
    
    print(f"Images saved as {unique_svg_filename} and {unique_jpg_filename}")
    
    # 返回文件名部分（中间的InChIKey或带下标的文件名）
    # 这里返回不带扩展名的文件名
    return os.path.splitext(unique_svg_filename)[0]

def process_smiles_file(file_path, smiles_column, output_dir, highlight_atoms, atom_colors):
    # 创建JPG和SVG输出目录
    svg_dir = os.path.join(output_dir, "svg")
    jpg_dir = os.path.join(output_dir, "jpg")
    if not os.path.exists(svg_dir):
        os.makedirs(svg_dir)
    if not os.path.exists(jpg_dir):
        os.makedirs(jpg_dir)
    
    # 读取SMILES文件
    if file_path.endswith('.csv'):
        df = pd.read_csv(file_path)
    elif file_path.endswith('.txt'):
        df = pd.read_csv(file_path, delimiter='\t')
    else:
        raise ValueError("Only CSV or TXT files are supported")
    
    # 存储生成的文件名
    file_names = []

    # 迭代处理每个SMILES
    for idx, smiles in enumerate(df[smiles_column].dropna()):
        short_filename = draw_and_save_smiles(smiles, svg_dir, jpg_dir, idx + 1, highlight_atoms, atom_colors)
        if short_filename:
            file_names.append(short_filename)
        else:
            file_names.append("Invalid_SMILES")
    
    # 将生成的文件名添加到CSV的最后一列
    df['Generated_Filename'] = file_names
    
    # 将修改后的CSV保存到原文件位置
    output_csv_path = os.path.join(output_dir, 'top_50_smiles_3p3h.csv')
    df.to_csv(output_csv_path, index=False)
    print(f"Updated CSV saved as {output_csv_path}")

# 示例使用
smiles_column = 'smiles'  # 你的文件中的SMILES列名
file_path = '/home/tim/hmt/del_picture/smiles/top_50_smiles_3p3h.txt'  # 你的CSV或TXT文件路径
output_dir = '/home/tim/hmt/del_picture/smiles'  # 输出文件夹

# 可以在这里快速设置需要高亮的原子和颜色
process_smiles_file(file_path, smiles_column, output_dir, HIGHLIGHT_ATOMS, ATOM_COLORS)
```