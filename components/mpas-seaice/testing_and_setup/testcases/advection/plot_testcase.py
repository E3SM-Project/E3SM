from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
import numpy as np

#-------------------------------------------------------------

def plot_subfigure(axis, array, nCells, nEdgesOnCell, verticesOnCell, xCell, yCell, zCell, xVertex, yVertex, zVertex, cmin, cmax, cmap):

    xMin =  1.0e30
    xMax = -1.0e30
    yMin =  1.0e30
    yMax = -1.0e30

    cmap = plt.get_cmap(cmap)

    patches = []
    colors = []

    for iCell in range(0,nCells):

        if (yCell[iCell] > 0.0):

            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append((xVertex[iVertex],zVertex[iVertex]))

            colors.append(array[iCell])

            patches.append(Polygon(vertices))

            xMin = min(xMin,xVertex[iVertex])
            xMax = max(xMax,xVertex[iVertex])

            yMin = min(yMin,zVertex[iVertex])
            yMax = max(yMax,zVertex[iVertex])

    pc = PatchCollection(patches, cmap=cmap)
    pc.set_array(np.array(colors))
    pc.set_clim(cmin, cmax)

    axis.add_collection(pc)

    axis.set_xlim(xMin,xMax)
    axis.set_ylim(yMin,yMax)
    axis.set_aspect("equal")

    axis.ticklabel_format(style='plain')
    axis.tick_params(axis='x', \
                     which='both', \
                     bottom=False, \
                     top=False, \
                     labelbottom=False)
    axis.tick_params(axis='y', \
                     which='both', \
                     left=False, \
                     right=False, \
                     labelleft=False)

#-------------------------------------------------------------

def plot_testcase():

    nGrids = [2562,10242,40962,163842]
    testTypes = ["cosine_bell","slotted_cylinder"]
    methods = ["IR","IR","upwind"]
    iTimes = [0,-1,-1]

    for nGrid in nGrids:

        print("nGrid: ", nGrid)

        fig, axes = plt.subplots(3,4)

        iTestType = -1
        for testType in testTypes:
            iTestType += 1

            print("  Test type: ", testType)

            iMethod = -1
            for method, iTime in zip(methods,iTimes):
                iMethod += 1

                print("    Method: ", method, ", iTime: ", iTime)

                filenamein = "./output_%s_%s_%i/output.2000.nc" %(method,testType,nGrid)

                filein = Dataset(filenamein,"r")

                nCells = len(filein.dimensions["nCells"])

                nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
                verticesOnCell = filein.variables["verticesOnCell"][:]
                xCell = filein.variables["xCell"][:]
                yCell = filein.variables["yCell"][:]
                zCell = filein.variables["zCell"][:]
                xVertex = filein.variables["xVertex"][:]
                yVertex = filein.variables["yVertex"][:]
                zVertex = filein.variables["zVertex"][:]
                verticesOnCell[:] = verticesOnCell[:] - 1

                iceAreaCategory = filein.variables["iceAreaCategory"][:]

                filein.close()

                iceAreaCell = np.sum(iceAreaCategory,axis=(2,3))

                plot_subfigure(axes[iMethod,iTestType*2], iceAreaCell[iTime], nCells, nEdgesOnCell, verticesOnCell, xCell, yCell, zCell, xVertex, yVertex, zVertex, 0.0, 1.0, "viridis")

                iceAreaCellDiff = iceAreaCell[iTime] - iceAreaCell[0]

                if (iMethod != 0):
                    plot_subfigure(axes[iMethod,iTestType*2+1], iceAreaCellDiff, nCells, nEdgesOnCell, verticesOnCell, xCell, yCell, zCell, xVertex, yVertex, zVertex, -1.0, 1.0, "bwr")
                else:
                    axes[iMethod,iTestType*2+1].axis('off')

        plt.savefig("advection_%6.6i.png" %(nGrid),dpi=300)
        plt.cla()
        plt.close(fig)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_testcase()
