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

def strain_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    resolutions = [2562,10242,40962,163842]

    methods = ["wachspress", "pwl", "weak", "weakwachs"]



    lineColours = ["black","black","grey","grey"]
    lineStyles  = ["-","--","-","--"]

    latitudeLimit = 20.0

    xMin = 6e-3
    xMax = 1e-1

    # quadratic scaling
    #scale = 3e-5 / math.pow(xMin,2)
    #scaleMin = math.pow(xMin,2) * scale
    #scaleMax = math.pow(xMax,2) * scale

    # linear scaling
    scale = 2.5e-3 / math.pow(xMin,1)
    scaleMin = math.pow(xMin,1) * scale
    scaleMax = math.pow(xMax,1) * scale


    plt.figure(figsize=(4,3))

    plt.loglog([xMin, xMax],[scaleMin,scaleMax],linestyle=':', color='k')

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
            y.append(normE11)

        plt.loglog(x,y, marker='o', color=lineColours[iPlot], ls=lineStyles[iPlot], markersize=5.0)

        iPlot = iPlot + 1

    #legendLabels = ["Quadratic scaling","Wachspress", "PWL", "Weak", "WeakWachs"]
    legendLabels = ["Linear scaling","Wachspress", "PWL", "Weak", "WeakWachs"]

    plt.legend(legendLabels, frameon=False, loc=2, fontsize=8, handlelength=4)

    ax = plt.gca()
    ax.set_xlabel("Grid resolution")
    ax.set_ylabel(r"$L_2$ error norm")
    ax.set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_scaling.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_scaling()
