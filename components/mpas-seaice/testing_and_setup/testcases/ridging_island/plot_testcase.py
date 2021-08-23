from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

#-------------------------------------------------------------------------------

def plot_testcase():

    iTime = -1

    # read in file
    filein = Dataset("./output_island/output.2000.nc","r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])

    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]

    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell -= 1

    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex -= 1

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    uVelocity = filein.variables["uVelocity"][iTime,:]
    vVelocity = filein.variables["vVelocity"][iTime,:]

    iceAreaCell = filein.variables["iceAreaCell"][iTime,:]
    iceVolumeCell = filein.variables["iceVolumeCell"][iTime,:]

    filein.close()

    xmin = np.amin(xCell)
    xmax = np.amax(xCell)
    ymin = np.amin(yCell)
    ymax = np.amax(yCell)

    # get patches list
    patchesCell = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex],yVertex[iVertex]))
        patchesCell.append(Polygon(vertices,True))

    patchesVertex = []
    deletedVertices = []
    for iVertex in range(0,nVertices):
        vertices = []
        useVertex = True
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            if (iCell != -1):
                vertices.append((xCell[iCell],yCell[iCell]))
            else:
                useVertex = False
        if (useVertex):
            patchesVertex.append(Polygon(vertices,True))
        else:
            deletedVertices.append(iVertex)
    deletedVertices = np.array(deletedVertices)

    uVelocity = np.delete(uVelocity, deletedVertices)
    vVelocity = np.delete(vVelocity, deletedVertices)

    # patch collections
    pcUVelocity = PatchCollection(patchesVertex, cmap=plt.get_cmap("jet"))
    pcUVelocity.set_array(uVelocity)

    pcVVelocity = PatchCollection(patchesVertex, cmap=plt.get_cmap("jet"))
    pcVVelocity.set_array(vVelocity)

    pcIceAreaCell = PatchCollection(patchesCell, cmap=plt.get_cmap("jet"))
    pcIceAreaCell.set_array(iceAreaCell)

    pcIceVolumeCell = PatchCollection(patchesCell, cmap=plt.get_cmap("jet"))
    pcIceVolumeCell.set_array(iceVolumeCell)

    # plot
    fig, axes = plt.subplots(2,2)

    axes[0,0].add_collection(pcIceAreaCell)
    divider = make_axes_locatable(axes[0,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcIceAreaCell, cax=cax)
    axes[0,0].set_xlim((xmin,xmax))
    axes[0,0].set_ylim((ymin,ymax))
    axes[0,0].set_aspect('equal')
    axes[0,0].ticklabel_format(style='plain')
    axes[0,0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    axes[0,0].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    axes[0,1].add_collection(pcIceVolumeCell)
    divider = make_axes_locatable(axes[0,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcIceVolumeCell, cax=cax)
    axes[0,1].set_xlim((xmin,xmax))
    axes[0,1].set_ylim((ymin,ymax))
    axes[0,1].set_aspect('equal')
    axes[0,1].ticklabel_format(style='plain')
    axes[0,1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    axes[0,1].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    axes[1,0].add_collection(pcUVelocity)
    divider = make_axes_locatable(axes[1,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcUVelocity, cax=cax)
    axes[1,0].set_xlim((xmin,xmax))
    axes[1,0].set_ylim((ymin,ymax))
    axes[1,0].set_aspect('equal')
    axes[1,0].ticklabel_format(style='plain')
    axes[1,0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    axes[1,0].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    axes[1,1].add_collection(pcVVelocity)
    divider = make_axes_locatable(axes[1,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcVVelocity, cax=cax)
    axes[1,1].set_xlim((xmin,xmax))
    axes[1,1].set_ylim((ymin,ymax))
    axes[1,1].set_aspect('equal')
    axes[1,1].ticklabel_format(style='plain')
    axes[1,1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    axes[1,1].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    plt.savefig("ridging_island.png",dpi=300)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_testcase()
