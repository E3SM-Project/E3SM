from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

#--------------------------------------------------------

def L2_norm_vertex(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > latitudeLimit * degreesToRadians):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_vertex(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    latVertex = fileMPAS.variables["latVertex"][:]

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    strain11varAvg = fileMPAS.variables["strain11varAvg"][0,:]
    strain22varAvg = fileMPAS.variables["strain22varAvg"][0,:]
    strain12varAvg = fileMPAS.variables["strain12varAvg"][0,:]

    normE11 = L2_norm_vertex(strain11varAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normE22 = L2_norm_vertex(strain22varAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normE12 = L2_norm_vertex(strain12varAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_cell(numerical, analytical, nCells, latCell, areaCell, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (math.fabs(latCell[iCell]) > latitudeLimit * degreesToRadians):

            norm  = norm  + areaCell[iCell] * math.pow(numerical[iCell] - analytical[iCell],2)

            denom = denom + areaCell[iCell] * math.pow(analytical[iCell],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_cell(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    latCell = fileMPAS.variables["latCell"][:]

    areaCell = fileMPAS.variables["areaCell"][:]

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    strain22weak = fileMPAS.variables["strain22weak"][0,:]
    strain12weak = fileMPAS.variables["strain12weak"][0,:]

    normE11 = L2_norm_cell(strain11weak, strain11CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    normE22 = L2_norm_cell(strain22weak, strain22CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    normE12 = L2_norm_cell(strain12weak, strain12CellAnalytical, nCells, latCell, areaCell, latitudeLimit)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def get_resolution(filename, latitudeLimit):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    degreesToRadians = math.pi / 180.0

    dcEdge = fileMPAS.variables["dcEdge"][:]
    latEdge = fileMPAS.variables["latEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        if (math.fabs(latEdge[iEdge]) > latitudeLimit * degreesToRadians):
            resolution = resolution + dcEdge[iEdge]
            denom = denom + 1.0

    resolution = resolution / denom

    fileMPAS.close()

    return resolution

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

    axis.loglog([xMin, xMax], [scaleMinLin,  scaleMaxLin],  linestyle=':', color='k')
    axis.loglog([xMin, xMax], [scaleMinQuad, scaleMaxQuad], linestyle=':', color='k')

#--------------------------------------------------------

def strain_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    strains = ["strain11","strain22","strain12"]

    resolutions = [2562,10242,40962,163842]

    methods = ["wachspress", "pwl", "weak"]

    lineColours = ["black","grey","red"]

    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    latitudeLimit = 20.0

    fig, axes = plt.subplots(3,1,figsize=(5,10))

    iStrain = 0
    for strain in strains:

        xMin = 4e-2
        xMax = 8e-2

        # linear scaling
        scale = 1e-3 / math.pow(xMin,1)
        scaleMinLin = math.pow(xMin,1) * scale
        scaleMaxLin = math.pow(xMax,1) * scale

        # quadratic scaling
        scale = 1e-3 / math.pow(xMin,2)
        scaleMinQuad = math.pow(xMin,2) * scale
        scaleMaxQuad = math.pow(xMax,2) * scale

        axes[iStrain].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
        axes[iStrain].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

        iPlot = 0
        for method in methods:

            x = []
            y = []

            for resolution in resolutions:

                filename = "./output_%s_%i/output.2000.nc" %(method,resolution)
                filenameIC = "./ic_%i.nc" %(resolution)

                print(filename, filenameIC)

                if (method == "weak"):
                    normE11, normE22, normE12 = get_norm_cell(filenameIC, filename, latitudeLimit)
                else:
                    normE11, normE22, normE12 = get_norm_vertex(filenameIC, filename, latitudeLimit)

                x.append(get_resolution(filename, latitudeLimit))
                if (strain == "strain11"):
                    y.append(normE11)
                elif (strain == "strain22"):
                    y.append(normE22)
                elif (strain == "strain12"):
                    y.append(normE12)

            axes[iStrain].loglog(x,y, marker='o', color=lineColours[iPlot], ls="solid", markersize=5.0, label="%s" %(method.capitalize()))

            iPlot = iPlot + 1

        axes[iStrain].legend(frameon=False, loc=2, fontsize=8, handlelength=4)

        axes[iStrain].set_xlabel("Grid resolution")
        axes[iStrain].set_ylabel(r"$L_2$ error norm")
        #axes[iStrain].set_xlim([xMin, xMax])

        iStrain = iStrain + 1


    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_scaling()
