---
title: pymol蛋白大小测量
published: 2025-01-15
description: "如何使用pymol显示蛋白/小分子的长宽高"
image: "./122680362_p0.png"
tags: [pymol,图]
category: "指南"
draft: false
lang: ""
---


:::tip
准备工作 蛋白质 `pdb`，配体 `mol2`，画图软件 `pymol`
:::


# 前言
本文介绍了如何使用 PyMOL 显示蛋白/小分子的长宽高

# 流程
```python 
from pymol import cmd, cgo
from pymol.cgo import *
import numpy
from random import randint

def matriz_inercia(selection):
	model = cmd.get_model(selection)
	totmass = 0.0
	x,y,z = 0,0,0
	for a in model.atom:
		m = a.get_mass()
		x += a.coord[0]*m
		y += a.coord[1]*m
		z += a.coord[2]*m
		totmass += m
	global cM
	cM = numpy.array([x/totmass, y/totmass, z/totmass])


	I = []
	for index in range(9):
		I.append(0)

	for a in model.atom:
		temp_x, temp_y, temp_z = a.coord[0], a.coord[1], a.coord[2]
		temp_x -= x
		temp_y -= y
		temp_z -= z

		I[0] += a.get_mass() * (temp_y**2 + temp_z**2)
		I[1] -= a.get_mass() * temp_x * temp_y
		I[2] -= a.get_mass() * temp_x * temp_z
		I[3] -= a.get_mass() * temp_x * temp_y
		I[4] += a.get_mass() * (temp_x**2 + temp_z**2)
		I[5] -= a.get_mass() * temp_y * temp_z
		I[6] -= a.get_mass() * temp_x * temp_z
		I[7] -= a.get_mass() * temp_y * temp_z
		I[8] += a.get_mass() * (temp_x**2 + temp_y**2)

	global tensor
	tensor = numpy.array([(I[0:3]), (I[3:6]), (I[6:9])])

	global autoval, autovect, ord_autoval, ord_autovect
	autoval, autovect = numpy.linalg.eig(tensor)
	auto_ord = numpy.argsort(autoval)
	ord_autoval = autoval[auto_ord]
	ord_autovect_complete = autovect[:, auto_ord].T
	ord_autovect = numpy.around(ord_autovect_complete, 3)

	return ord_autoval

def translacion_cM(selection):
	model = cmd.get_model(selection)
	totmass = 0.0
	x,y,z = 0,0,0
	for a in model.atom:
		m = a.get_mass()
		x += a.coord[0]*m
		y += a.coord[1]*m
		z += a.coord[2]*m
		totmass += m
	cM = numpy.array([x/totmass, y/totmass, z/totmass])
	trans_array = ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -cM[0], -cM[1], -cM[2], 1])
	model_trans = cmd.transform_selection(selection, trans_array)


def rotacion_orig(selection):
	translacion_cM(selection)
	matriz_inercia(selection)
	global transf, transf_array, ord_autovect_array, transf_array_print
	ord_autovect_array = numpy.array([[ord_autovect[0][0], ord_autovect[0][1], ord_autovect[0][2]],
	                                  [ord_autovect[1][0], ord_autovect[1][1], ord_autovect[1][2]],
					  [ord_autovect[2][0], ord_autovect[2][1], ord_autovect[2][2]]])
	if numpy.linalg.det(ord_autovect_array) < 0.:
		ord_autovect_array = numpy.array([[ord_autovect[2][0], ord_autovect[2][1], ord_autovect[2][2]],
	                                  	  [ord_autovect[1][0], ord_autovect[1][1], ord_autovect[1][2]],
	 				  	  [ord_autovect[0][0], ord_autovect[0][1], ord_autovect[0][2]]])
	transf = numpy.transpose(ord_autovect_array)
	transf_array = numpy.array([transf[0][0], transf[0][1], transf[0][2], 0,
	                            transf[1][0], transf[1][1], transf[1][2], 0,
				    transf[2][0], transf[2][1], transf[2][2], 0,
				    0, 0, 0, 1])


def transformar(selection):
	rotacion_orig(selection)
	model_rot = cmd.transform_selection(selection, transf_array, homogenous=0, transpose=1);

def draw_IABB(selection):
	transformar(selection)

	([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

	print("The Inertia Axis Aligned Bounding Box (IABB) dimensions are (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))
	print("The Inertia Axis Aligned Bounding Box (IABB) volume is %.2f A3" % ((maxX-minX)*(maxY-minY)*(maxZ-minZ)))


	minX = minX
	minY = minY
	minZ = minZ
	maxX = maxX
	maxY = maxY
	maxZ = maxZ

	boundingBox = [
                LINEWIDTH, float(2),

                BEGIN, LINES,
                COLOR, float(1), float(0), float(0),

                VERTEX, minX, minY, minZ,
                VERTEX, minX, minY, maxZ,
                VERTEX, minX, maxY, minZ,
                VERTEX, minX, maxY, maxZ,
                VERTEX, maxX, minY, minZ,
                VERTEX, maxX, minY, maxZ,
                VERTEX, maxX, maxY, minZ,
                VERTEX, maxX, maxY, maxZ,
                VERTEX, minX, minY, minZ,
                VERTEX, maxX, minY, minZ,
                VERTEX, minX, maxY, minZ,
                VERTEX, maxX, maxY, minZ,
                VERTEX, minX, maxY, maxZ,
                VERTEX, maxX, maxY, maxZ,
                VERTEX, minX, minY, maxZ,
                VERTEX, maxX, minY, maxZ,
                VERTEX, minX, minY, minZ,
                VERTEX, minX, maxY, minZ,
                VERTEX, maxX, minY, minZ,
                VERTEX, maxX, maxY, minZ,
                VERTEX, minX, minY, maxZ,
                VERTEX, minX, maxY, maxZ,
                VERTEX, maxX, minY, maxZ,
                VERTEX, maxX, maxY, maxZ,

                END
        ]

	p4 = '_4' + str(randint(0, 100))
	p5 = '_5' + str(randint(0, 100))
	p6 = '_6' + str(randint(0, 100))
	p7 = '_7' + str(randint(0, 100))
	cmd.pseudoatom (pos=[minX, minY, minZ], object=p4)
	cmd.pseudoatom (pos=[minX, minY, maxZ], object=p5)
	cmd.pseudoatom (pos=[minX, maxY, minZ], object=p6)
	cmd.pseudoatom (pos=[maxX, minY, minZ], object=p7)
	cmd.distance(None, p4, p7)
	cmd.distance(None, p4, p6)
	cmd.distance(None, p4, p5)
	cmd.hide("nonbonded")

	boxName = "box_IABB_" + str(randint(0, 100))
	cmd.load_cgo(boundingBox,boxName)
	return boxName

def draw_3DSIZE(selection):
	draw_IABB(selection)

cmd.extend ("draw_3DSIZE", draw_3DSIZE)

```

## 运行
命令行操作
```bash
run draw_Protein.py
draw_3DSIZE object     # object为Pymol中蛋白的名字
draw_3DSIZE sele       
```
