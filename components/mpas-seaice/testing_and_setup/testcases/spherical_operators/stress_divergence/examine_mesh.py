from netCDF4 import Dataset
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Plot vertices and edges of a cell')

parser.add_argument('-i', dest='filenameIn', required=True, help='Input mesh file')
parser.add_argument('-c', dest='iCell', required=True, type=int, help='iCell to plot')

args = parser.parse_args()


filein = Dataset(args.filenameIn,"r")

nCells = len(filein.dimensions["nCells"])
nVertices = len(filein.dimensions["nVertices"])
nEdges = len(filein.dimensions["nEdges"])

xCell = filein.variables["xCell"][:]
yCell = filein.variables["yCell"][:]
zCell = filein.variables["zCell"][:]

xVertex = filein.variables["xVertex"][:]
yVertex = filein.variables["yVertex"][:]
zVertex = filein.variables["zVertex"][:]

xEdge = filein.variables["xEdge"][:]
yEdge = filein.variables["yEdge"][:]
zEdge = filein.variables["zEdge"][:]

nEdgesOnCell = filein.variables["nEdgesOnCell"][:]

edgesOnCell = filein.variables["edgesOnCell"][:]
verticesOnCell = filein.variables["verticesOnCell"][:]

edgesOnCell[:] -= 1
verticesOnCell[:] -= 1

#iCell = 123
iCell = args.iCell

colors = ["red","green","blue","orange","yellow","purple","magenta","cyan","teal"]

fig = plt.figure()
ax = plt.axes(projection="3d")

for iEdgeOnCell in range(0,nEdgesOnCell[iCell]):

    iVertexOnCell1 = iEdgeOnCell
    iVertexOnCell2 = iEdgeOnCell + 1
    if (iVertexOnCell2 == nEdgesOnCell[iCell]): iVertexOnCell2 = 0

    iVertex1 = verticesOnCell[iCell,iVertexOnCell1]
    iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

    iEdge = edgesOnCell[iCell,iEdgeOnCell]

    print(xVertex[iVertex1],yVertex[iVertex1],zVertex[iVertex1])
    print(xEdge  [iEdge],   yEdge  [iEdge],   zEdge  [iEdge])
    print(xVertex[iVertex2],yVertex[iVertex2],zVertex[iVertex2])
    print()

    ax.scatter3D([xVertex[iVertex1]],[yVertex[iVertex1]],[zVertex[iVertex1]], c=colors[iEdgeOnCell], marker="+");
    ax.scatter3D([xEdge  [iEdge]],   [yEdge  [iEdge]],   [zEdge  [iEdge]],    c=colors[iEdgeOnCell], marker=".");
    ax.scatter3D([xVertex[iVertex2]],[yVertex[iVertex2]],[zVertex[iVertex2]], c=colors[iEdgeOnCell], marker="x");

plt.show()

filein.close()
