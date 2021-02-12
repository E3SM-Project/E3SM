from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#--------------------------------------------------------

def L2_norm_weak(numerical, analytical, nCells, areaCell, useCell):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (useCell[iCell]):

            norm  = norm  + areaCell[iCell] * math.pow(numerical[iCell] - analytical[iCell],2)

            denom = denom + areaCell[iCell] * math.pow(analytical[iCell],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_weak(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    areaCell = fileMPAS.variables["areaCell"][:]

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    strain22weak = fileMPAS.variables["strain22weak"][0,:]
    strain12weak = fileMPAS.variables["strain12weak"][0,:]

    useCell = get_use_cell(filename)

    normE11 = L2_norm_weak(strain11weak, strain11CellAnalytical, nCells, areaCell, useCell)
    normE22 = L2_norm_weak(strain22weak, strain22CellAnalytical, nCells, areaCell, useCell)
    normE12 = L2_norm_weak(strain12weak, strain12CellAnalytical, nCells, areaCell, useCell)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_var(numerical, analytical, nVertices, areaTriangle, useVertex):

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

def get_norm_var(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    latVertex = fileMPAS.variables["latVertex"][:]

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    strain11Vertex = fileMPAS.variables["strain11varAvg"][0,:]
    strain22Vertex = fileMPAS.variables["strain22varAvg"][0,:]
    strain12Vertex = fileMPAS.variables["strain12varAvg"][0,:]

    useVertex = get_use_vertex(filename)

    normE11 = L2_norm_var(strain11Vertex, strain11VertexAnalytical, nVertices, areaTriangle, useVertex)
    normE22 = L2_norm_var(strain22Vertex, strain22VertexAnalytical, nVertices, areaTriangle, useVertex)
    normE12 = L2_norm_var(strain12Vertex, strain12VertexAnalytical, nVertices, areaTriangle, useVertex)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def get_resolution(filename):

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

def get_use_cell(filenameIn):

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

    useCell = np.ones(nCells,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            useCell[iCell] = 0
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                useCell[iCell2] = 0

    return useCell

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

def strain_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    strains = ["strain11","strain22","strain12"]

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

    fig, axes = plt.subplots(3,1,figsize=(5,10))

    iStrain = 0
    for strain in strains:

        xMin = 6e-3
        xMax = 1e-2

        # linear scaling
        scale = 4e-4 / math.pow(xMin,1)
        scaleMinLin = math.pow(xMin,1) * scale
        scaleMaxLin = math.pow(xMax,1) * scale

        # quadratic scaling
        scale = 4e-4 / math.pow(xMin,2)
        scaleMinQuad = math.pow(xMin,2) * scale
        scaleMaxQuad = math.pow(xMax,2) * scale

        axes[iStrain].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
        axes[iStrain].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

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

                    if (operatorMethod == "wachspress" or
                        operatorMethod == "pwl"):
                        normE11, normE22, normE12 = get_norm_var(filenameIC, filename)
                    else:
                        normE11, normE22, normE12 = get_norm_weak(filenameIC, filename)

                    x.append(get_resolution(filename))
                    if (strain == "strain11"):
                        y.append(normE11)
                    elif (strain == "strain22"):
                        y.append(normE22)
                    elif (strain == "strain12"):
                        y.append(normE12)

                axes[iStrain].loglog(x, y, marker='o', color=lineColours[iPlot], ls=lineStyles[gridType], markersize=5.0, label="%s" %(operatorMethod.capitalize()))

                iPlot = iPlot + 1

        axes[iStrain].legend(frameon=False, loc=2, fontsize=8, handlelength=4)

        axes[iStrain].set_xlabel("Grid resolution")
        axes[iStrain].set_ylabel(r"$L_2$ error norm")
        #axes[iStrain].set_xlim([xMin, xMax])
        axes[iStrain].set_title("%s" %(strain))

        iStrain = iStrain + 1


    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_scaling()
