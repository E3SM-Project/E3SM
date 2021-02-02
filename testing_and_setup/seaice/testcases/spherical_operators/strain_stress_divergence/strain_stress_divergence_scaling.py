from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

#--------------------------------------------------------

def L2_norm(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

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

def get_norm(filenameIC, filename, latitudeLimit):

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

    normU = L2_norm(stressDivergenceU, stressDivergenceUAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normV = L2_norm(stressDivergenceV, stressDivergenceVAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)

    fileMPAS.close()

    return normU, normV

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

def strain_stress_divergence_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    resolutions = [2562,10242,40962,163842]
    methods = ["wachspress","pwl","weak","weakwachs","weakpwl","wachspress_alt","pwl_alt","weakwachs_alt","weakpwl_alt"]

    lineColours = ["black","black","grey","red","red","blue","blue","green","green"]
    lineStyles  = ["-","--","-","-","--","-","--","-","--"]

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

            normU, normV = get_norm(filenameIC, filename, latitudeLimit)

            x.append(get_resolution(filename, latitudeLimit))
            y.append(normU)

        plt.loglog(x,y, marker='o', color=lineColours[iPlot], ls=lineStyles[iPlot], markersize=5.0)

        iPlot = iPlot + 1

    #legendLabels = ["Quadratic scaling","Wachspress", "PWL", "Weak", "WeakWachs"]
    legendLabels = ["Linear scaling","Wachspress", "PWL", "Weak", "WeakWachs", "WeakPWL", "Wachspress Alt", "PWL Alt","WeakWachs Alt","WeakPWL Alt"]

    plt.legend(legendLabels, frameon=False, loc=2, fontsize=8, handlelength=4)

    ax = plt.gca()
    ax.set_xlabel("Grid resolution")
    ax.set_ylabel(r"$L_2$ error norm")
    ax.set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_stress_divergence_scaling.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_stress_divergence_scaling()
