---
title: 突出苯磺酰胺的smiles分子图绘制
published: 2024-09-27
description: "如何根据几个选定的颜色绘制还算好看的smlies图，并突出苯磺酰胺基团"
image: "./222.jpg"
tags: [python,分子,图]
category: "画图"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 分子 `smiles`

:::

## 图片绘制代码
```python 
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import os
from rdkit.Chem import inchi

# 控制开关
HIGHLIGHT_O_F = True  # 设置为False可以不高亮O和F原子

# 颜色设置（可以根据需要修改）
ATOM_COLORS = {
    'N': (1.0, 0.7569, 0.6235),  # FFC19F
    'S': (0.9725, 0.7059, 0.8627),  # F8B4DC
    'O': (0.0000, 0.7882, 0.7176) ,   # #008578 深青绿色
    'F': (0.9725, 0.7059, 0.8627), # F8B4DC
}
BENZENESULFONAMIDE_COLOR = (0.9922, 0.5765, 0.7686)  # #fd93c4 亮粉红色
HIGHLIGHT_BOND_COLOR =  (1.0000, 0.7941, 0.9725)  # #ffe4f8 非常浅的粉红色

# SMILES string of the molecule
smiles = "COC(=O)[C@H](CCC(=O)N[C@@H](CC(=O)O)C(=O)OCc1ccccc1)NC(=O)c1ccc(S(N)(=O)=O)cc1"

# Create RDKit molecule object from SMILES
mol = Chem.MolFromSmiles(smiles)

# Generate 2D coordinates for the molecule
AllChem.Compute2DCoords(mol)

def assign_colors_by_atom_type(mol):
    color_map = {}
    for i, atom in enumerate(mol.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        if atom_symbol in ATOM_COLORS:
            if atom_symbol in ['O', 'F'] and not HIGHLIGHT_O_F:
                continue
            color_map[i] = ATOM_COLORS[atom_symbol]
    return color_map

# Function to identify benzenesulfonamide group
def identify_benzenesulfonamide(mol):
    pattern = Chem.MolFromSmarts('c1ccc(S(=O)(N)=O)cc1')
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return list(matches[0])
    return []

# Drawing settings
draw_options = rdMolDraw2D.MolDrawOptions()
draw_options.useBWAtomPalette()  # Turn off default coloring
draw_options.fixedBondLength = 30

# Prepare the drawer and color map
drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
drawer.SetDrawOptions(draw_options)
color_map = assign_colors_by_atom_type(mol)

# Identify benzenesulfonamide group
benzenesulfonamide_atoms = identify_benzenesulfonamide(mol)

# Add benzenesulfonamide atoms to the color map
for atom in benzenesulfonamide_atoms:
    color_map[atom] = BENZENESULFONAMIDE_COLOR

# Combine N, S atoms and benzenesulfonamide atoms, and O, F if HIGHLIGHT_O_F is True
highlight_atoms = list(set([i for i, atom in enumerate(mol.GetAtoms()) 
                            if atom.GetSymbol() in ['N', 'S'] or 
                            (HIGHLIGHT_O_F and atom.GetSymbol() in ['O', 'F'])] + 
                           benzenesulfonamide_atoms))

# Identify bonds in the benzenesulfonamide group for bolding
highlight_bonds = []
for i, bond in enumerate(mol.GetBonds()):
    if bond.GetBeginAtomIdx() in benzenesulfonamide_atoms and bond.GetEndAtomIdx() in benzenesulfonamide_atoms:
        highlight_bonds.append(i)

# Draw the molecule
drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightAtomColors=color_map,
                    highlightBonds=highlight_bonds, highlightBondColors={i: HIGHLIGHT_BOND_COLOR for i in highlight_bonds})

drawer.FinishDrawing()
svg_output = drawer.GetDrawingText()

# Generate InChIKey for the molecule
inchi_key = inchi.MolToInchiKey(mol)

# Create a base filename using the InChIKey
base_filename = f"molecule_{inchi_key[:8]}"  # Using first 8 characters of InChIKey

# Save as SVG
svg_filename = f"{base_filename}.svg"
with open(svg_filename, "w") as svg_file:
    svg_file.write(svg_output.replace('svg:', ''))

# Save as JPG
jpg_filename = f"{base_filename}.jpg"
img = Draw.MolToImage(mol, size=(900, 900), highlightAtoms=highlight_atoms, highlightAtomColors=color_map,
                      highlightBonds=highlight_bonds)
img.save(jpg_filename)

print(f"Molecule SMILES: {smiles}")
print(f"Molecule Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"Images saved as {svg_filename} and {jpg_filename}")
```