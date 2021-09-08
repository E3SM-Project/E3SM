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

def stress_divergence_scaling():

    resolutions = [2562,10242,40962,163842]

    methods = ["wachspress", "pwl", "weak", "wachspress_alt", "pwl_alt"]


    latitudeLimit = 20.0

    # plot options
    legendLabels = {"wachspress":"Wachs.",
                    "pwl":"PWL",
                    "weak":"Weak",
                    "wachspress_alt":"Wachs. alt",
                    "pwl_alt":"PWL alt."}

    lineColors = {"wachspress":"black",
                  "pwl":"grey",
                  "weak":"red",
                  "wachspress_alt":"black",
                  "pwl_alt":"grey"}

    lineStyles = {"wachspress":"solid",
                  "pwl":"solid",
                  "weak":"solid",
                  "wachspress_alt":"dashed",
                  "pwl_alt":"dashed"}

    markers = {"wachspress":"+",
               "pwl":"x",
               "weak":"^",
               "wachspress_alt":"+",
               "pwl_alt":"x"}

    # plot
    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    # scaling lines
    xMin = 2.5e-2
    xMax = 5e-2
    yMin = 2.5e-3

    fig, axes = plt.subplots(figsize=(4,3))

    scaling_lines(axes, xMin, xMax, yMin)

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

        axes.loglog(x,y, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])

        iPlot = iPlot + 1

    axes.legend(frameon=False, loc=2, fontsize=8, handlelength=4)

    axes.set_xlabel("Grid resolution")
    axes.set_ylabel(r"$L_2$ error norm")
    #axes.set_xlim([6e-3, 1e-1])
    axes.set_title(r'(b) $(\nabla \cdot \sigma)_u$')
    axes.set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes.set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes.set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes.set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_divergence_scaling.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_divergence_scaling()
