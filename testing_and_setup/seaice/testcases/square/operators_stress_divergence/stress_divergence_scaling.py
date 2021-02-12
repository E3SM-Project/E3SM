from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#--------------------------------------------------------

def L2_norm(numerical, analytical, nVertices, areaTriangle, useVertex):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (useVertex[iVertex]):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm(filenameIC, filename, useVertex):

    fileIC = Dataset(filenameIC, "r")

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    latVertex = fileMPAS.variables["latVertex"][:]

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    stressDivergenceU = fileMPAS.variables["stressDivergenceU"][0,:]
    stressDivergenceV = fileMPAS.variables["stressDivergenceV"][0,:]

    normU = L2_norm(stressDivergenceU, stressDivergenceUAnalytical, nVertices, areaTriangle, useVertex)
    normV = L2_norm(stressDivergenceV, stressDivergenceVAnalytical, nVertices, areaTriangle, useVertex)

    fileMPAS.close()

    return normU, normV

#--------------------------------------------------------

def get_resolution(filename, useVertex):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    degreesToRadians = math.pi / 180.0

    dcEdge = fileMPAS.variables["dcEdge"][:]
    latEdge = fileMPAS.variables["latEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        resolution = resolution + dcEdge[iEdge]
        denom = denom + 1.0

    resolution = resolution / denom

    fileMPAS.close()

    return resolution

#--------------------------------------------------------

def get_use_vertex(filenameIn):

    fileIn = Dataset(filenameIn,"r")
    nCells = len(fileIn.dimensions["nCells"])
    nVertices = len(fileIn.dimensions["nVertices"])
    nEdgesOnCell = fileIn.variables["nEdgesOnCell"][:]
    cellsOnCell = fileIn.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1
    verticesOnCell = fileIn.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1
    interiorCell = fileIn.variables["interiorCell"][0,:]
    fileIn.close()

    useVertex = np.ones(nVertices,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                if (iCell2 < nCells):
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell2]):
                        iVertex = verticesOnCell[iCell2,iVertexOnCell]
                        useVertex[iVertex] = 0

    return useVertex

#--------------------------------------------------------

def stress_divergence_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    operatorMethods = ["wachspress","pwl","weak"]

    gridTypes = ["hex","quad"]
    #gridTypes = ["quad"]

    grids = {"hex" :["0082x0094",
                     "0164x0188",
                     "0328x0376",
                     "0656x0752"],
             "quad":["0080x0080",
                     "0160x0160",
                     "0320x0320",
                     "0640x0640"]}
    #grids = {"hex" :["0082x0094"],
    #         "quad":["0080x0080"]}


    lineColours = ["black","grey","red"]

    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    #markers = {"hex":"o",
    #           "quad":"s"}
    markers = ["P","X","*"]


    xMin = 6e-3
    xMax = 1e-2

    # quadratic scaling
    #scale = 3e-5 / math.pow(xMin,2)
    #scaleMin = math.pow(xMin,2) * scale
    #scaleMax = math.pow(xMax,2) * scale

    # linear scaling
    scale = 2.5e-3 / math.pow(xMin,1)
    scaleMin = math.pow(xMin,1) * scale
    scaleMax = math.pow(xMax,1) * scale


    plt.figure(figsize=(4,3))

    #plt.loglog([xMin, xMax],[scaleMin,scaleMax],linestyle=':', color='k')


    for gridType in gridTypes:

        print("Grid type: ", gridType)

        iPlot = 0
        for operatorMethod in operatorMethods:

            print("  Operator Method: ", operatorMethod)

            x = []
            y = []

            for grid in grids[gridType]:

                print("    Grid: ", grid)

                filenameIC = "ic_%s_%s.nc" %(gridType,grid)
                filename = "output_%s_%s_%s/output.2000.nc" %(gridType, operatorMethod, grid)

                print("      ", filename, filenameIC)

                useVertex = get_use_vertex(filename)

                normU, normV = get_norm(filenameIC, filename, useVertex)

                x.append(get_resolution(filename, useVertex))
                y.append(normU)

            print(x,y)
            plt.loglog(x, y, color=lineColours[iPlot], ls=lineStyles[gridType], marker=markers[iPlot], markersize=5.0)

            iPlot = iPlot + 1

    #legendLabels = ["Quadratic scaling","Wachspress", "PWL", "Weak", "WeakWachs"]
    legendLabels = ["Wachspress", "PWL", "Weak"]

    plt.legend(legendLabels, frameon=False, loc=2, fontsize=8, handlelength=4)

    ax = plt.gca()
    ax.set_xlabel("Grid resolution")
    ax.set_ylabel(r"$L_2$ error norm")
    #ax.set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_divergence_scaling.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_divergence_scaling()
