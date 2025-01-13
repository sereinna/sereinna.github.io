---
title: pymol蛋白图自动绘制
published: 2025-01-13
description: "如何使用pymol自动绘制好看的蛋白质-配体相互作用图"
image: "./125226470_p0.png"
tags: [pymol,图]
category: "指南"
draft: false
lang: ""
---


:::tip
准备工作 蛋白质 `pdb`，配体 `mol2`，画图软件 `pymol`
:::


# 前言
本文介绍了如何使用 PyMOL 绘制蛋白质的远景图和近景图，包括设置显示模式、调整颜色和参数，以及标记相互作用残基等步骤。

# 流程
```python 
#!/usr/env python

from pymol import cmd 
from pymol.cmd import util
from pymol import stored
import pymol2
import os,sys
import time
import subprocess



#hydrogen bonds extend
def list_hb(selection,selection2=None,cutoff=3,angle=55,mode=1,hb_list_name='hbonds'):
    """
    USAGE

    list_hb selection, [selection2 (default=None)], [cutoff (default=3.2)],
                       [angle (default=55)], [mode (default=1)],
                       [hb_list_name (default='hbonds')]

    The script automatically adds a requirement that atoms in the
    selection (and selection2 if used) must be either of the elements N or
    O.

    If mode is set to 0 instead of the default value 1, then no angle
    cutoff is used, otherwise the angle cutoff is used and defaults to 55
    degrees.

    e.g.
    To get a list of all H-bonds within chain A of an object
      list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds

    To get a list of H-bonds between chain B and everything else:
      list_hb 1tl9 & c. b, 1tl9 &! c. b

    """
    cutoff=float(cutoff)
    angle=float(angle)
    mode=float(mode)
    # ensure only N and O atoms are in the selection
    selection3 = selection + " & e. h & ((neighbor e. n) |(neighbor e. o))"
    if not selection2:
        selection4 = selection + " & e. n+o"
    else:
        selection4 = selection2 + " & e. n+o"
    hb1 = cmd.find_pairs(selection3,selection4,mode=mode,cutoff=cutoff,angle=angle)

    selection3 = selection + " & e. n+o" 
    if not selection2:
        selection4 = selection + " & e. h & (neighbor e. n+o)"
    else:
        selection4 = selection2 + " & e. h & (neighbor e. n+o)"
    hb2=cmd.find_pairs(selection3,selection4,mode=mode,cutoff=cutoff,angle=angle)
    
    for pairs in hb1+hb2:
        #cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'print("%1s/%3s`%s/%-4s " % (chain,resn,resi,name))')
        #cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'print("%1s/%3s`%s/%-4s " % (chain,resn,resi,name))')
        cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
        stored.list=[]
        cmd.iterate(str(pairs[0][0])+' and index '+str(pairs[0][1]),"stored.list.append(name)")
     
        hb_atom=stored.list[0]
        if hb_atom == 'O' or  hb_atom == 'HN' or  hb_atom == 'H':
            
            cmd.select('bb_'+hb_list_name,'(byres ('+pairs[0][0]+' and index '+str(pairs[0][1])+')) and ((name CA) or (name C) or (name O) or (name N) or (name NH) or (name H))')
            cmd.show("sticks",'bb_'+hb_list_name)  
            cmd.set("stick_radius",0.1,'bb_'+hb_list_name)

#Halogen bonds extend
def list_xb(selection,selection2=None,cutoff=3.5,angle=55,mode=1,xb_list_name='xbonds'):
    """
    USAGE

    list_xb selection, [selection2 (default=None)], [cutoff (default=3.2)],
                       [angle (default=55)], [mode (default=1)],
                       [hb_list_name (default='xbonds')]

    The script automatically adds a requirement that atoms in the
    selection (and selection2 if used) must be either of the elements N or
    O.

    If mode is set to 0 instead of the default value 1, then no angle
    cutoff is used, otherwise the angle cutoff is used and defaults to 55
    degrees.


    """ 
    cutoff=float(cutoff)
    angle=float(angle)
    mode=float(mode)
    # ensure only N and O atoms are in the selection
    selection = selection + " & e. n+o+s"
    if not selection2:
        selection3 = selection + " & e. cl+br+i"
    else:
        selection3 = selection2 + " & e. cl+br+i"
    xb = cmd.find_pairs(selection,selection3,mode=mode,cutoff=cutoff,angle=angle)

    
    # sort the list for easier reading
    #hb.sort(lambda x,y:(cmp(x[0][1],y[0][1])))

    for pairs in xb:
        #cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'print("%1s/%3s`%s/%-4s " % (chain,resn,resi,name))')
        #cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'print("%1s/%3s`%s/%-4s " % (chain,resn,resi,name))')
        cmd.distance(xb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
        stored.list=[]
        cmd.iterate(str(pairs[0][0])+' and index '+str(pairs[0][1]),"stored.list.append(name)")
      
        xb_atom=stored.list[0]
        if xb_atom == 'O' or  xb_atom == 'N':
            
            cmd.select('bb_'+xb_list_name,'(byres ('+pairs[0][0]+' and index '+str(pairs[0][1])+')) and ((name CA) or (name C) or (name O) or (name N) or (name NH) or (name H) )')
            cmd.show("sticks",'bb_'+xb_list_name)   
            cmd.set("stick_radius",0.1,'bb_'+xb_list_name)


#pi interations bonds extend
def pi_interations(ligand='ligand',selection='ligand',selection2='site',pp_list_name='ppbonds'):
 
    from rdkit import Chem

    stored.list=[]
    cmd.iterate(selection2,"stored.list.append(resn+' '+resi+' '+chain)")
    resd_names=set(stored.list)
    mol = Chem.MolFromMolFile(ligand)
    ssr = Chem.GetSymmSSSR(mol)
    for ring in ssr:
        j=0
       
        #b=" ".join(str(x) for x in list(ring))
        #print(b)
        for i in ring:
            atom=mol.GetAtomWithIdx(i)
            if atom.GetIsAromatic() : j += 1 
          
        if j == len(ring):
                k=1
                for resd_name_str in resd_names:
                 
                    resd_name=resd_name_str.split(' ')
                   
                    if  resd_name[0] == 'PHE':
                        pp_distance=cmd.distance('pipi_PHE_'+str(i)+'_'+str(k),selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+') and ((resn PHE) and (name CG+CD1+CD2+CE1+CD2+CZ))',cutoff=5,mode=4 )
                        if pp_distance > 4.2:
                            cmd.delete('pipi_PHE_'+str(i)+'_'+str(k))    
                 #    cmd.distance('pipi_PHE_'+i+'_'+j,'ligand and index '+'+'.join(ring), 'site and ((resn. PHE) and (name CG+CD1+CD2+CE1+CD2+CZ)) and byres (ligand around 4)',"5","4" )
                    elif resd_name[0] == 'HIS' or resd_name[0] == 'HID' or resd_name[0] == 'HIE':
                       pp_distance=cmd.distance('pipi_HIS_'+str(i)+'_'+str(k),selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+') and ((resn HIS+HIE+HID) and (name CG+ND1+CD2+CE1+CE2))',cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_HIS_'+str(i)+'_'+str(k))
                   #    cmd.distance('pipi_HIS_'+i+'_'+j,'ligand and index '+list(ring).split(',','+'), '((resn. HIS+HIE+HID) and (name CG+ND1+CD2+CE1+CE2)) and byres (ligand around 4)',"5","4" )
                    elif resd_name[0] == 'TYR':                        
                       pp_distance=cmd.distance('pipi_TYR_'+str(i)+'_'+str(k),selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+')  and ((resn TYR)and(name CG+CD1+CD2+CE1+CD2+CZ))',cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_TYR_'+str(i)+'_'+str(k))
                    elif resd_name[0] == 'TRP' :
                    
                       pp_distance=cmd.distance('pipi_TRP1_'+str(i)+'_'+str(k),selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+')  and ((resn TRP) and (name CG+CD1+CD2+NE1+CE2))',cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_TRP1_'+str(i)+'_'+str(k))
                       pp_distance=cmd.distance('pipi_TRP2_'+str(i)+'_'+str(k),selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+')  and ((resn TRP) and (name CD2+NE1+CE2+CE3+CZ2+CZ3))',cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_TPR2_'+str(i)+'_'+str(k))
                   #    cmd.distance('pipi_TRP1_'+i+'_'+j,'ligand and index '+list(ring).split(',','+'), '(resn. TRP) and (name CG+CD1+CD2+NE1+CE2) ) and byres (ligand around 4)',"5","4" )
                   #    cmd.distance('pipi_TRP2_'+i+'_'+j,'ligand and index '+list(ring).split(',','+'), '((resn. TRP) and (name CD2+NE1+CE2+CE3+CZ2+CZ3)) and byres (ligand around 4)',"5","4" )
                    elif resd_name[0] == 'LYS' :
                        
                        #cmd.iterate('(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+' )  and ((resn LYS) and (name NZ))','resi+resn')
                        cation_pi_distance=cmd.distance('cation_pi_LYS_'+str(i)+'_'+str(k),selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+' ) and (resn LYS) and (name NZ)',cutoff=5,mode=4 )                        
                        if cation_pi_distance > 5:
                            cmd.delete('cation_pi_LYS_'+str(i)+'_'+str(k))
                          

                    elif resd_name[0] == 'ARG' :
                      
                        #cmd.iterate('(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+' )  and ((resn LYS) and (name CZ))','resi+resn')
                        cation_pi_distance=cmd.distance('cation_pi_ARG_'+str(i)+'_'+str(k),selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+' ) and (resn LYS) and (name CZ)',cutoff=5,mode=4 )                        
                        if cation_pi_distance > 5:
                            cmd.delete('cation_pi_ARG_'+str(i)+'_'+str(k))
                            
                    k += 1
 
    stored.list=[]
    cmd.iterate(selection+' and (formal_charge >0) and (elem N) ',"stored.list.append(index)")
    charge_idxs=set(stored.list)
    for charge_idx in charge_idxs:
        
        #b=" ".join(str(x) for x in list(ring))
        #print(b)
        k=1
        for resd_name_str in resd_names:
   
            resd_name=resd_name_str.split(' ')
    
            if  resd_name[0] == 'PHE':
                cation_pi_distance=cmd.distance('cation_pi_PHE_'+str(charge_idx)+'_'+str(k),selection1=selection+' and index '+str(charge_idx),selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+') and ((resn PHE) and (name CG+CD1+CD2+CE1+CD2+CZ))',cutoff=5,mode=4 )
                if cation_pi_distance > 5:
                    cmd.delete('cation_pi_PHE_'+str(charge_idx)+'_'+str(k))    
         #    cmd.distance('pipi_PHE_'+i+'_'+j,'ligand and index '+'+'.join(ring), 'site and ((resn. PHE) and (name CG+CD1+CD2+CE1+CD2+CZ)) and byres (ligand around 4)',"5","4" )
            elif resd_name[0] == 'HIS' or resd_name[0] == 'HID' or resd_name[0] == 'HIE':
               cation_pi_distance=cmd.distance('cation_pi_HIS_'+str(charge_idx)+'_'+str(k),selection1=selection+' and index '+str(charge_idx),selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+') and ((resn HIS+HIE+HID) and (name CG+ND1+CD2+CE1+CE2))',cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_HIS_'+str(charge_idx)+'_'+str(k))
           #    cmd.distance('pipi_HIS_'+i+'_'+j,'ligand and index '+list(ring).split(',','+'), '((resn. HIS+HIE+HID) and (name CG+ND1+CD2+CE1+CE2)) and byres (ligand around 4)',"5","4" )
            elif resd_name[0] == 'TYR':                        
               cation_pi_distance=cmd.distance('cation_pi_TYR_'+str(charge_idx)+'_'+str(k),selection1=selection+' and index '+str(charge_idx),selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+')  and ((resn TYR)and(name CG+CD1+CD2+CE1+CD2+CZ))',cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_TYR_'+str(charge_idx)+'_'+str(k))
            elif resd_name[0] == 'TRP' :
            
               cation_pi_distance=cmd.distance('cation_pi_TRP1_'+str(charge_idx)+'_'+str(k),selection1=selection+' and index '+str(charge_idx),selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+')  and ((resn TRP) and (name CG+CD1+CD2+NE1+CE2))',cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_TRP1_'+str(charge_idx)+'_'+str(k))
               cation_pi_distance=cmd.distance('cation_pi_TRP2_'+str(charge_idx)+'_'+str(k),selection1=selection+' and index '+str(charge_idx),selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+')  and ((resn TRP) and (name CD2+NE1+CE2+CE3+CZ2+CZ3))',cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_TPR2_'+str(charge_idx)+'_'+str(k))
           #    cmd.distance('pipi_TRP1_'+i+'_'+j,'ligand and index '+list(ring).split(',','+'), '(resn. TRP) and (name CG+CD1+CD2+NE1+CE2) ) and byres (ligand around 4)',"5","4" )
           #    cmd.distance('pipi_TRP2_'+i+'_'+j,'ligand and index '+list(ring).split(',','+'), '((resn. TRP) and (name CD2+NE1+CE2+CE3+CZ2+CZ3)) and byres (ligand around 4)',"5","4" )
            k += 1

    cmd.group(pp_list_name,'open')
    cmd.group(pp_list_name,'pipi_*')
    cmd.group(pp_list_name,'cation_pi_*')
    cmd.group(pp_list_name,'close')                      

#---main---

#---protrin---

def cpi(cpxname):
    command="set_name %s, pdb"%(cpxname)
    cmd.do(command)
    name="lig"
    cmd.do('extract ligand,organic')
    #___cattoon___
    cmd.show("cartoon",'pdb')
    cmd.set("cartoon_transparency",0.5)
    cmd.set("cartoon_oval_width","0.15")
    cmd.set("cartoon_oval_length","0.8")
    cmd.set("sphere_scale",0.3)


    #___site___
    cmd.group('g_site','open')
    cmd.select("site","(byres (ligand around 4)) and (not name N) and (not name O) and (not name C) ")
    cmd.select("site_ca","(name CA) and (byres (ligand around 4))")
    #cmd.select("ligand2","(byres (ligand around 4))")
    #cmd.select("ligand3","(name CA) and ligand2")
#    cmd.set("cartoon_flat_cycles",0)
    cmd.select("protein_c","(ele c) and pdb")
    #cmd.select("ligand_c","(ele c) and ligand")
    cmd.select("site_c","(ele c) and site")
    #cmd.color("tv_green","protein_c")
    cmd.color("gray80","protein_c")
    #cmd.color("cyan","ligand_c")
#    cmd.color("tv_green","site_c")
    cmd.label("site_ca","resn+resi")
    #cmd.set("valence","0")
    #cmd.show("lines","site")
    #cmd.set("line_width",4)
    cmd.show("sticks","site")
    cmd.set("stick_radius",0.1,'site')
    cmd.group('g_site','site')
    cmd.group('g_site','site_ca')
    cmd.group('g_site','protein_c')
    cmd.group('g_site','site_c')
    cmd.group('g_site','close')

    #___surface___
    cmd.group("g_surface","open")
    util.protein_vacuum_esp('pdb',mode=2,quiet=0,_self=cmd)
    cmd.set("surface_carve_selection","pdb")
    cmd.set("surface_carve_cutoff",6)
#    cmd.set("transparency",0,"pdb_e_chg")
#    cmd.set("surface_solvent",1)
    cmd.group("g_surface","pdb_e_map")
    cmd.group("g_surface","pdb_e_chg")
    cmd.group("g_surface","pdb_e_pot")
    cmd.group("g_surface","close")


    #---first ligand---
    cmd.group("g1_"+name,"open")
    cmd.show("sticks","ligand")
    #cmd.distance("hbonds",selection1="ligand",selection2="pdb",cutoff=3.5,mode=2)
    #___hysrogen bond___
    list_hb(selection="pdb",selection2="ligand",hb_list_name="hbonds")
    #___Pi interactions___
    #pi_interations(ligand="ligand",selection='ligand',selection2='site',pp_list_name='ppbonds')
    #___Salt bridge positive___
    cmd.distance('sbp', "ligand and  (formal_charge >0) and (elem N)", "pdb and (resn ASP+Glu) and (name OD*+OE*)", "5.5", "0")
    #___Salt bridge NEGATIVE___
    cmd.distance('sbn', "ligand and  (formal_charge < 0) and (elem O)", "pdb and  (resn Lys and name NZ) or (resn arg and name NE+NH*)", "5.5", "0")
    #___Halogen bond__
    list_xb(selection="pdb",selection2="ligand",xb_list_name='xbonds')

    #cmd.group("g1_"+name,"ligand")
    cmd.group("interaction_"+name,"hbonds")
    #cmd.group("interaction_"+name,"pi")
    cmd.group("interaction_"+name,"sbp")
    cmd.group("interaction_"+name,"sbn")
    cmd.group("interaction_"+name,"ppbonds")
    cmd.group("interaction_"+name,"bb_hbonds")
    cmd.group("interaction_"+name,"bb_xbonds")
    cmd.group("interaction_"+name,"xbonds")
    cmd.group("interaction_"+name,"close")
            
    #---general---   
    cmd.hide('sticks',"h. and (ele c extend 1)")
    cmd.set("label_color","black")
    cmd.set("label_size","16")
    cmd.set ("label_font_id","5")
    cmd.set("stick_h_scale","1")

    #cmd.hide("lines","ele h")

    #cmd.hide("cartoon")

    cmd.color_deep("yellow", 'hbond*', 0)
    cmd.color_deep("deeppurple", 'xbond*', 0)
    cmd.color_deep("aquamarine", 'ppbonds*', 0)
    cmd.color_deep("violet", 'sbp*', 0)
    cmd.color_deep("violet", 'sbn*', 0)
    cmd.bg_color("white")
    cmd.disable("g_surface")
    cmd.enable("pdb")
    cmd.do('set ray_trace_mode, 1')
    cmd.do('set ray_shadows,0')
    cmd.do('set specular, 0')
    cmd.do('space cmyk')
    cmd.do('set ray_trace_color, [0,0,0]')
    cmd.do('set stick_h_scale,1')

cmd.extend("cpi", cpi)


```

## 结构主体显示
### pdb
- **功能**：蛋白质结构整体
- **内容**：包含主链和侧链的完整三维结构
- **显示**：默认以`cartoon`方式显示
- **说明**：可调整透明度和颜色

### ligand
- **功能**：配体分子整体
- **内容**：从复合物中提取的有机小分子
- **显示**：默认以`stick`方式显示
- **说明**：可独立调整显示方式

## 结合位点分析 (g_site)

### site
- **功能**：配体4Å范围内的结合口袋
- **内容**：包含所有可能的相互作用残基
- **显示**：以`stick`方式显示侧链
- **范围**：`ligand around 4`

### site_ca
- **功能**：结合位点残基的α碳
- **显示**：显示残基标签(序号和类型)
- **用途**：标识关键残基位置
- **标签**：`resn+resi`

### protein_c
- **功能**：结合位点区域的蛋白质碳原子
- **用途**：确定蛋白质骨架位置
- **颜色**：设置特定颜色以区分

### site_c
- **功能**：结合位点残基的碳原子
- **用途**：分析疏水相互作用
- **显示**：突出显示关键接触

## 分子表面显示 (g_surface)

### pdb_e_chg
- **功能**：静电势着色的分子表面
- **显示**：电荷分布可视化
- **颜色**：
  * 红色：负电势
  * 蓝色：正电势

### pdb_e_map
- **功能**：分子表面网格
- **用途**：显示表面轮廓
- **说明**：用于可视化表面特征

### pdb_e_pot
- **功能**：静电势等势面
- **显示**：电势分布梯度
- **内容**：包含电势值标度

## 相互作用分析 (interaction_lig)

### hbonds
- **功能**：氢键相互作用
- **显示**：供体-受体对
- **参数**：默认距离阈值`3.5Å`
- **颜色**：黄色

### bb_hbonds
- **功能**：主链氢键
- **说明**：与蛋白质主链形成的氢键
- **重要性**：对结合稳定性重要

### sbp (Salt Bridge Positive)
- **功能**：正电荷盐桥
- **类型**：配体正电荷与蛋白质负电荷残基
- **参数**：距离阈值`5.5Å`
- **颜色**：紫色

### sbn (Salt Bridge Negative)
- **功能**：负电荷盐桥
- **类型**：配体负电荷与蛋白质正电荷残基
- **参数**：距离阈值`5.5Å`
- **颜色**：紫色


