---
title: smiles分子图绘制
published: 2024-09-18
description: "如何根据几个选定的颜色绘制还算好看的smlies图"
image: "./121711656_p0.png"
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

# SMILES string of the molecule
smiles = "CC(C)(C)C(=O)Nc1ccc(S(=O)(=O)N)cc1"

# Create RDKit molecule object from SMILES
mol = Chem.MolFromSmiles(smiles)

# Generate 2D coordinates for the molecule
AllChem.Compute2DCoords(mol)

# Define colors by atom type in hexadecimal converted to RGB
atom_colors = {
    'N': (1.0, 0.7569, 0.6235),  # FFC19F
    'O': (1.0, 0.7098, 0.7686),  # FFB5C4
    'F': (0.9725, 0.7059, 0.8627), # F8B4DC
    'P': (1.0, 0.8549, 0.4863),  # FFDA7C
}

def assign_colors_by_atom_type(mol):
    color_map = {}
    for i, atom in enumerate(mol.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        if atom_symbol in atom_colors:
            color_map[i] = atom_colors[atom_symbol]
    return color_map

# Drawing settings
draw_options = rdMolDraw2D.MolDrawOptions()
draw_options.useBWAtomPalette()  # Turn off default coloring
draw_options.fixedBondLength = 30

# Prepare the drawer and color map
drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
drawer.SetDrawOptions(draw_options)
color_map = assign_colors_by_atom_type(mol)
non_carbon_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() != 'C']

# Draw the molecule
if non_carbon_atoms:
    drawer.DrawMolecule(mol, highlightAtoms=non_carbon_atoms, highlightAtomColors=color_map)
else:
    drawer.DrawMolecule(mol)

drawer.FinishDrawing()
svg_output = drawer.GetDrawingText()

# Save as SVG
with open("molecule11.svg", "w") as svg_file:
    svg_file.write(svg_output.replace('svg:', ''))

# Save as JPG
img = Draw.MolToImage(mol, size=(300, 300), highlightAtoms=non_carbon_atoms, highlightColorDict=color_map)
img.save("molecule11.jpg")

print(f"Molecule SMILES: {smiles}")
print(f"Molecule Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print("Images saved as molecule.svg and molecule.jpg")
```