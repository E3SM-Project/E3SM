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

    cellsOnEdge = fileMPAS.variables["cellsOnEdge"][:]
    cellsOnEdge[:] = cellsOnEdge[:] - 1

    dcEdge = fileMPAS.variables["dcEdge"][:]
    latEdge = fileMPAS.variables["latEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        iCell1 = cellsOnEdge[iEdge,0]
        iCell2 = cellsOnEdge[iEdge,1]
        if (iCell1 != -1 and iCell2 != -1):
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

def scaling_lines(axis, xMin, xMax, yMin):

    # linear scaling
    scale = yMin / math.pow(xMin,1)
    scaleMinLin = math.pow(xMin,1) * scale
    scaleMaxLin = math.pow(xMax,1) * scale

    # quadratic scaling
    scale = yMin / math.pow(xMin,2)
    scaleMinQuad = math.pow(xMin,2) * scale
    scaleMaxQuad = math.pow(xMax,2) * scale

    axis.loglog([xMin, xMax], [scaleMinLin,  scaleMaxLin],  linestyle=':', color='k', label='_nolegend_')
    axis.loglog([xMin, xMax], [scaleMinQuad, scaleMaxQuad], linestyle=':', color='k', label='_nolegend_')

#--------------------------------------------------------

def strain_stress_divergence_scaling():

    # options
    operatorMethods = ["wachspress","pwl","weak","weakwachs","weakpwl","wachsavg","pwlavg"]

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

    stressDivergences = ["U","V"]


    # plot options
    lineColours = {"wachspress":"black",
                   "pwl":"grey",
                   "weak":"red",
                   "weakwachs":"blue",
                   "weakpwl":"green",
                   "wachsavg":"purple",
                   "pwlavg":"orange"}

    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    legendLabels = {"wachspress":"Wachs.",
                    "pwl":"PWL",
                    "weak":"Weak",
                    "weakwachs":"Weak Wachs.",
                    "weakpwl":"Weak PWL",
                    "wachsavg":"Wachs. Avg.",
                    "pwlavg":"PWL Avg."}

    stressDivergenceLabels = {"U":r"(a) $(\nabla \cdot \sigma)_u$",
                              "V":r"(b) $(\nabla \cdot \sigma)_v$"}


    # scaling lines
    xMin = 2e-3
    xMax = 3.5e-3
    yMin = 9e-3

    # plot
    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(2, 1, figsize=(5,7))

    j = 0
    for stressDivergence in stressDivergences:

        scaling_lines(axes[j], xMin, xMax, yMin)

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
                    if (stressDivergence == "U"):
                        y.append(normU)
                    elif (stressDivergence == "V"):
                        y.append(normV)

                if (gridType == "hex"):
                    legendLabel = legendLabels[operatorMethod]
                else:
                    legendLabel = "_nolegend_"
                axes[j].loglog(x, y, marker='o', color=lineColours[operatorMethod], ls=lineStyles[gridType], markersize=5.0, label=legendLabel)

                iPlot = iPlot + 1

        axes[j].legend(frameon=False, loc=4, fontsize=8, handlelength=4)

        axes[j].set_xlabel("Grid resolution")
        axes[j].set_ylabel(r"$L_2$ error norm")
        axes[j].set_title(stressDivergenceLabels[stressDivergence])

        j += 1

    plt.tight_layout()
    plt.savefig("strain_stress_divergence_scaling.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_stress_divergence_scaling()
