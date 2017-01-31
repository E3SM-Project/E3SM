#!/usr/bin/env python
import os
import numpy as np

import argparse

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--file", dest="filename", help="Path to grid file", metavar="FILE", required=True)
parser.add_argument("-w", "--weights", dest="weight_field", help="Field to weight block partition file on.", metavar="VAR")

args = parser.parse_args()

if not args.weight_field:
	print "Weight field missing. Defaulting to unweighted graphs."
	weighted_parts = False
else:
	weighted_parts = True

dev_null = open(os.devnull, 'w')
grid = NetCDFFile(args.filename, 'r')

nCells = len(grid.dimensions['nCells'])
nEdges = len(grid.dimensions['nEdges'])

nEdgesOnCell = grid.variables['nEdgesOnCell'][:]
cellsOnCell = grid.variables['cellsOnCell'][:] - 1
if weighted_parts:
	try:
		weights = grid.variables[args.weight_field][:]
	except:
		print args.weight_field, ' not found in file. Defaulting to un-weighted partitions.'
		weighted_parts = False
grid.close()

nEdges = 0
for i in np.arange(0, nCells):
	for j in np.arange(0,nEdgesOnCell[i]):
		if cellsOnCell[i][j] != -1:
			nEdges = nEdges + 1

nEdges = nEdges/2

graph = open('graph.info', 'w+')
if weighted_parts:
	graph.write('%s %s 010\n'%(nCells, nEdges))
else:
	graph.write('%s %s\n'%(nCells, nEdges))

for i in np.arange(0, nCells):
	if weighted_parts:
		graph.write('%s '%int(weights[i]))

	for j in np.arange(0,nEdgesOnCell[i]):
		if(cellsOnCell[i][j] >= 0):
			graph.write('%s '%(cellsOnCell[i][j]+1))
	graph.write('\n')
graph.close()

