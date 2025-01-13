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


#卤键扩展函数
#通过设定距离和角度阈值来寻找可能形成卤键的原子对,并将它们可视化显示。对于含N/O的残基,还会显示其主链骨架结构,便于观察卤键在蛋白质结构中的位置。
def list_xb(selection,selection2=None,cutoff=3.5,angle=55,mode=1,xb_list_name='xbonds'):
    """
    使用方法

    list_xb selection, [selection2 (默认=None)], [cutoff (默认=3.2)],
                       [angle (默认=55)], [mode (默认=1)],
                       [xb_list_name (默认='xbonds')]

    此脚本自动要求选区(selection和selection2如果使用)中的原子必须是
    N或O元素。

    如果mode设置为0(默认值为1),则不使用角度阈值限制,否则使用默认
    55度的角度阈值。
    """ 
    cutoff=float(cutoff)
    angle=float(angle)
    mode=float(mode)
    # 确保选区中只包含N和O原子
    selection = selection + " & e. n+o+s"
    if not selection2:
        selection3 = selection + " & e. cl+br+i"
    else:
        selection3 = selection2 + " & e. cl+br+i"
    xb = cmd.find_pairs(selection,selection3,mode=mode,cutoff=cutoff,angle=angle)
    
    # 对列表排序以便阅读
    #hb.sort(lambda x,y:(cmp(x[0][1],y[0][1])))

    for pairs in xb:
        # 创建卤键distance对象
        cmd.distance(xb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
        stored.list=[]
        cmd.iterate(str(pairs[0][0])+' and index '+str(pairs[0][1]),"stored.list.append(name)")
      
        xb_atom=stored.list[0]
        if xb_atom == 'O' or  xb_atom == 'N':
            # 选择并显示包含N/O原子的残基主链
            cmd.select('bb_'+xb_list_name,'(byres ('+pairs[0][0]+' and index '+str(pairs[0][1])+')) and ((name CA) or (name C) or (name O) or (name N) or (name NH) or (name H) )')
            cmd.show("sticks",'bb_'+xb_list_name)   
            cmd.set("stick_radius",0.1,'bb_'+xb_list_name)




#π相互作用键扩展函数
#用于识别和显示蛋白质-配体复合物中的π-π堆积和阳离子-π相互作用。首先使用RDKit识别配体中的芳香环系统,然后检查这些环与蛋白质中芳香氨基酸(PHE、TYR、TRP、HIS)以及带正电荷的氨基酸(LYS、ARG)之间可能存在的相互作用，代码使用距离阈值(4.2Å用于π-π堆积,5Å用于阳离子-π相互作用)来筛选有效的相互作用。最后将所有识别到的相互作用归类并可视化显示。
def pi_interations(ligand='ligand',selection='ligand',selection2='site',pp_list_name='ppbonds'):
    """
    功能：识别并显示分子间的π-π堆积和阳离子-π相互作用
    
    参数说明：
    ligand: 配体的mol文件路径
    selection: 第一个选择对象(通常为配体)
    selection2: 第二个选择对象(通常为结合位点)
    pp_list_name: 相互作用列表名称
    
    注意：
    - π-π堆积使用4.2Å作为距离阈值
    - 阳离子-π相互作用使用5Å作为距离阈值
    """
 
    from rdkit import Chem

    # 获取结合位点中的残基信息
    stored.list=[]
    cmd.iterate(selection2,"stored.list.append(resn+' '+resi+' '+chain)")
    resd_names=set(stored.list)
    
    # 使用RDKit读取配体并识别环系统
    mol = Chem.MolFromMolFile(ligand)
    ssr = Chem.GetSymmSSSR(mol)
    
    # 遍历所有环系统
    for ring in ssr:
        j=0
       
        # 检查环的芳香性
        for i in ring:
            atom=mol.GetAtomWithIdx(i)
            if atom.GetIsAromatic() : j += 1 
          
        # 如果是完全芳香环，检查可能的相互作用
        if j == len(ring):
                k=1
                for resd_name_str in resd_names:
                    resd_name=resd_name_str.split(' ')
                   
                    # 处理PHE残基的π-π相互作用
                    if  resd_name[0] == 'PHE':
                        pp_distance=cmd.distance('pipi_PHE_'+str(i)+'_'+str(k),
                            selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),
                            selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+') and ((resn PHE) and (name CG+CD1+CD2+CE1+CD2+CZ))',
                            cutoff=5,mode=4 )
                        if pp_distance > 4.2:
                            cmd.delete('pipi_PHE_'+str(i)+'_'+str(k))    
                 
                    # 处理HIS残基的π-π相互作用
                    elif resd_name[0] == 'HIS' or resd_name[0] == 'HID' or resd_name[0] == 'HIE':
                       pp_distance=cmd.distance('pipi_HIS_'+str(i)+'_'+str(k),
                           selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),
                           selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+') and ((resn HIS+HIE+HID) and (name CG+ND1+CD2+CE1+CE2))',
                           cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_HIS_'+str(i)+'_'+str(k))
                   
                    # 处理TYR残基的π-π相互作用
                    elif resd_name[0] == 'TYR':                        
                       pp_distance=cmd.distance('pipi_TYR_'+str(i)+'_'+str(k),
                           selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),
                           selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+')  and ((resn TYR)and(name CG+CD1+CD2+CE1+CD2+CZ))',
                           cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_TYR_'+str(i)+'_'+str(k))
                    
                    # 处理TRP残基的π-π相互作用（两个环）
                    elif resd_name[0] == 'TRP' :
                       # 五元环
                       pp_distance=cmd.distance('pipi_TRP1_'+str(i)+'_'+str(k),
                           selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),
                           selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+')  and ((resn TRP) and (name CG+CD1+CD2+NE1+CE2))',
                           cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_TRP1_'+str(i)+'_'+str(k))
                       # 六元环
                       pp_distance=cmd.distance('pipi_TRP2_'+str(i)+'_'+str(k),
                           selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),
                           selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+')  and ((resn TRP) and (name CD2+NE1+CE2+CE3+CZ2+CZ3))',
                           cutoff=5,mode=4 )
                       if pp_distance > 4.2:
                           cmd.delete('pipi_TPR2_'+str(i)+'_'+str(k))
                   
                    # 处理LYS残基的阳离子-π相互作用
                    elif resd_name[0] == 'LYS' :
                        cation_pi_distance=cmd.distance('cation_pi_LYS_'+str(i)+'_'+str(k),
                            selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),
                            selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+' ) and (resn LYS) and (name NZ)',
                            cutoff=5,mode=4 )                        
                        if cation_pi_distance > 5:
                            cmd.delete('cation_pi_LYS_'+str(i)+'_'+str(k))
                          
                    # 处理ARG残基的阳离子-π相互作用
                    elif resd_name[0] == 'ARG' :
                        cation_pi_distance=cmd.distance('cation_pi_ARG_'+str(i)+'_'+str(k),
                            selection1=selection+' and index '+("+".join(str(x) for x in list(ring))),
                            selection2='(chain '+resd_name[2]+') and ( resi '+ resd_name[1]+' ) and (resn LYS) and (name CZ)',
                            cutoff=5,mode=4 )                        
                        if cation_pi_distance > 5:
                            cmd.delete('cation_pi_ARG_'+str(i)+'_'+str(k))
                            
                    k += 1
 
    # 处理带正电荷的氮原子的阳离子-π相互作用
    stored.list=[]
    cmd.iterate(selection+' and (formal_charge >0) and (elem N) ',"stored.list.append(index)")
    charge_idxs=set(stored.list)
    
    # 遍历所有带正电荷的氮原子
    for charge_idx in charge_idxs:
        k=1
        for resd_name_str in resd_names:
            resd_name=resd_name_str.split(' ')
    
            # 检查与PHE残基的相互作用
            if  resd_name[0] == 'PHE':
                cation_pi_distance=cmd.distance('cation_pi_PHE_'+str(charge_idx)+'_'+str(k),
                    selection1=selection+' and index '+str(charge_idx),
                    selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+') and ((resn PHE) and (name CG+CD1+CD2+CE1+CD2+CZ))',
                    cutoff=5,mode=4 )
                if cation_pi_distance > 5:
                    cmd.delete('cation_pi_PHE_'+str(charge_idx)+'_'+str(k))    
            
            # 检查与HIS残基的相互作用
            elif resd_name[0] == 'HIS' or resd_name[0] == 'HID' or resd_name[0] == 'HIE':
               cation_pi_distance=cmd.distance('cation_pi_HIS_'+str(charge_idx)+'_'+str(k),
                   selection1=selection+' and index '+str(charge_idx),
                   selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+') and ((resn HIS+HIE+HID) and (name CG+ND1+CD2+CE1+CE2))',
                   cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_HIS_'+str(charge_idx)+'_'+str(k))
           
            # 检查与TYR残基的相互作用
            elif resd_name[0] == 'TYR':                        
               cation_pi_distance=cmd.distance('cation_pi_TYR_'+str(charge_idx)+'_'+str(k),
                   selection1=selection+' and index '+str(charge_idx),
                   selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+')  and ((resn TYR)and(name CG+CD1+CD2+CE1+CD2+CZ))',
                   cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_TYR_'+str(charge_idx)+'_'+str(k))
            
            # 检查与TRP残基的相互作用
            elif resd_name[0] == 'TRP' :
               # 五元环
               cation_pi_distance=cmd.distance('cation_pi_TRP1_'+str(charge_idx)+'_'+str(k),
                   selection1=selection+' and index '+str(charge_idx),
                   selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+')  and ((resn TRP) and (name CG+CD1+CD2+NE1+CE2))',
                   cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_TRP1_'+str(charge_idx)+'_'+str(k))
               # 六元环
               cation_pi_distance=cmd.distance('cation_pi_TRP2_'+str(charge_idx)+'_'+str(k),
                   selection1=selection+' and index '+str(charge_idx),
                   selection2='(chain '+resd_name[2]+') and( resi '+ resd_name[1]+')  and ((resn TRP) and (name CD2+NE1+CE2+CE3+CZ2+CZ3))',
                   cutoff=5,mode=4 )
               if cation_pi_distance > 5:
                   cmd.delete('cation_pi_TPR2_'+str(charge_idx)+'_'+str(k))
            k += 1

    # 将所有相互作用归类到指定组
    cmd.group(pp_list_name,'open')
    cmd.group(pp_list_name,'pipi_*')
    cmd.group(pp_list_name,'cation_pi_*')
    cmd.group(pp_list_name,'close')                 

#---主程序---

#---蛋白质处理---

def cpi(cpxname):
    """
    复合物结构分析和可视化主函数
    
    参数：
    cpxname: 复合物结构文件名
    """
    # 设置结构名称
    command="set_name %s, pdb"%(cpxname)
    cmd.do(command)
    name="lig"
    # 提取有机小分子作为配体
    cmd.do('extract ligand,organic')
    
    #---卡通表示设置---
    cmd.show("cartoon",'pdb')
    cmd.set("cartoon_transparency",0.5)    # 设置卡通透明度
    cmd.set("cartoon_oval_width","0.15")   # 设置卡通椭圆宽度
    cmd.set("cartoon_oval_length","0.8")   # 设置卡通椭圆长度
    cmd.set("sphere_scale",0.3)            # 设置球体大小

    #---结合位点处理---
    cmd.group('g_site','open')
    # 选择配体4埃范围内的残基，排除主链原子
    cmd.select("site","(byres (ligand around 4)) and (not name N) and (not name O) and (not name C) ")
    # 选择结合位点的CA原子
    cmd.select("site_ca","(name CA) and (byres (ligand around 4))")
    # 选择碳原子
    cmd.select("protein_c","(ele c) and pdb")
    cmd.select("site_c","(ele c) and site")
    # 设置颜色
    cmd.color("gray80","protein_c")
    # 添加残基标签
    cmd.label("site_ca","resn+resi")
    # 显示结合位点棍状模型
    cmd.show("sticks","site")
    cmd.set("stick_radius",0.1,'site')
    # 组织结合位点显示
    cmd.group('g_site','site')
    cmd.group('g_site','site_ca')
    cmd.group('g_site','protein_c')
    cmd.group('g_site','site_c')
    cmd.group('g_site','close')

    #---表面处理---
    cmd.group("g_surface","open")
    # 计算蛋白质真空静电势
    util.protein_vacuum_esp('pdb',mode=2,quiet=0,_self=cmd)
    # 设置表面雕刻
    cmd.set("surface_carve_selection","pdb")
    cmd.set("surface_carve_cutoff",6)
    # 组织表面显示
    cmd.group("g_surface","pdb_e_map")
    cmd.group("g_surface","pdb_e_chg")
    cmd.group("g_surface","pdb_e_pot")
    cmd.group("g_surface","close")

    #---配体处理---
    cmd.group("g1_"+name,"open")
    cmd.show("sticks","ligand")
    
    #---相互作用分析---
    # 氢键分析
    list_hb(selection="pdb",selection2="ligand",hb_list_name="hbonds")
    # 正电荷盐桥
    cmd.distance('sbp', "ligand and  (formal_charge >0) and (elem N)", 
                "pdb and (resn ASP+Glu) and (name OD*+OE*)", "5.5", "0")
    # 负电荷盐桥
    cmd.distance('sbn', "ligand and  (formal_charge < 0) and (elem O)", 
                "pdb and  (resn Lys and name NZ) or (resn arg and name NE+NH*)", "5.5", "0")
    # 卤键分析
    list_xb(selection="pdb",selection2="ligand",xb_list_name='xbonds')

    #---相互作用分组---
    cmd.group("interaction_"+name,"hbonds")
    cmd.group("interaction_"+name,"sbp")
    cmd.group("interaction_"+name,"sbn")
    cmd.group("interaction_"+name,"ppbonds")
    cmd.group("interaction_"+name,"bb_hbonds")
    cmd.group("interaction_"+name,"bb_xbonds")
    cmd.group("interaction_"+name,"xbonds")
    cmd.group("interaction_"+name,"close")
            
    #---通用显示设置---   
    # 隐藏氢原子
    cmd.hide('sticks',"h. and (ele c extend 1)")
    # 标签设置
    cmd.set("label_color","black")
    cmd.set("label_size","16")
    cmd.set("label_font_id","5")
    # 棍状模型设置
    cmd.set("stick_h_scale","1")

    #---相互作用颜色设置---
    cmd.color_deep("yellow", 'hbond*', 0)      # 氢键显示为黄色
    cmd.color_deep("deeppurple", 'xbond*', 0)  # 卤键显示为深紫色
    cmd.color_deep("aquamarine", 'ppbonds*', 0) # π相互作用显示为海蓝色
    cmd.color_deep("violet", 'sbp*', 0)        # 盐桥显示为紫色
    cmd.color_deep("violet", 'sbn*', 0)        # 盐桥显示为紫色
    
    #---渲染设置---
    cmd.bg_color("white")                      # 设置白色背景
    cmd.disable("g_surface")                   # 默认关闭表面显示
    cmd.enable("pdb")                          # 启用蛋白质显示
    cmd.do('set ray_trace_mode, 1')           # 设置光线追踪模式
    cmd.do('set ray_shadows,0')               # 关闭阴影
    cmd.do('set specular, 0')                 # 关闭镜面反射
    cmd.do('space cmyk')                      # 使用CMYK色彩空间
    cmd.do('set ray_trace_color, [0,0,0]')    # 设置光线追踪颜色
    cmd.do('set stick_h_scale,1')             # 设置氢原子棍状模型比例

# 扩展PyMOL命令
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


