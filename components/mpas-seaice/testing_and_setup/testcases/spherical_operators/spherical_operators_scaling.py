from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#--------------------------------------------------------
# strain
#--------------------------------------------------------

def L2_norm_strain_vertex(numerical, analytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit):

    # variational

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            for iCellOnVertex in range(0,vertexDegree):
                iCell2 = cellsOnVertex[iVertex,iCellOnVertex]
                if (iCell == iCell2):
                    area = kiteAreasOnVertex[iVertex,iCellOnVertex]

            if (math.fabs(latVertex[iVertex]) > math.radians(latitudeLimit)):

                norm  = norm  + area * math.pow(numerical[iCell,iVertexOnCell] - analytical[iCell,iVertexOnCell],2)

                denom = denom + area * math.pow(analytical[iCell,iVertexOnCell],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_vertex(numerical, analytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit):

    # variational
    norm = 0.0

    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]

            if (math.fabs(latVertex[iVertex]) > math.radians(latitudeLimit)):
                norm = max(norm, math.fabs(numerical[iCell,iVertexOnCell] - analytical[iCell,iVertexOnCell]))

    return norm

#--------------------------------------------------------

def get_norm_strain_vertex(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # variational

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    kiteAreasOnVertex = fileGrid.variables["kiteAreasOnVertex"][:]

    verticesOnCell[:] = verticesOnCell[:] - 1
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    fileGrid.close()

    # IC
    fileIC = Dataset(filenameIC, "r")

    # nVertices
    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    #strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    #strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    strain11VarAnalytical = np.zeros((nCells, maxEdges))
    #strain22VarAnalytical = np.zeros((nCells, maxEdges))
    #strain12VarAnalytical = np.zeros((nCells, maxEdges))
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            strain11VarAnalytical[iCell,iVertexOnCell] = strain11VertexAnalytical[iVertex]
            #strain22VarAnalytical[iCell,iVertexOnCell] = strain22VertexAnalytical[iVertex]
            #strain12VarAnalytical[iCell,iVertexOnCell] = strain12VertexAnalytical[iVertex]

    # results
    fileMPAS = Dataset(filename, "r")

    # nCells, maxEdges
    strain11var = fileMPAS.variables["strain11var"][0,:]
    #strain22var = fileMPAS.variables["strain22var"][0,:]
    #strain12var = fileMPAS.variables["strain12var"][0,:]

    fileMPAS.close()

    if (normType == "l2"):
        normE11 = L2_norm_strain_vertex(strain11var, strain11VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_vertex(strain22var, strain22VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_vertex(strain12var, strain12VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_vertex(strain11var, strain11VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_vertex(strain22var, strain22VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_vertex(strain12var, strain12VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_strain_cell(numerical, analytical, nCells, latCell, areaCell, latitudeLimit):

    # weak

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (math.fabs(latCell[iCell]) > np.radians(latitudeLimit)):

            norm  = norm  + areaCell[iCell] * math.pow(numerical[iCell] - analytical[iCell],2)

            denom = denom + areaCell[iCell] * math.pow(analytical[iCell],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_cell(numerical, analytical, nCells, latCell, areaCell, latitudeLimit):

    # weak
    norm  = 0.0
    for iCell in range(0,nCells):

        if (math.fabs(latCell[iCell]) > np.radians(latitudeLimit)):
            norm = max(norm, math.fabs(numerical[iCell] - analytical[iCell]))

    return norm

#--------------------------------------------------------

def get_norm_strain_cell(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # weak

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    latCell = fileGrid.variables["latCell"][:]
    areaCell = fileGrid.variables["areaCell"][:]

    fileGrid.close()

    # IC
    fileIC = Dataset(filenameIC, "r")

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    #strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    #strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    #strain22weak = fileMPAS.variables["strain22weak"][0,:]
    #strain12weak = fileMPAS.variables["strain12weak"][0,:]

    if (normType == "l2"):
        normE11 = L2_norm_strain_cell(strain11weak, strain11CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_cell(strain22weak, strain22CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_cell(strain12weak, strain12CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_cell(strain11weak, strain11CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_cell(strain22weak, strain22CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_cell(strain12weak, strain12CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_strain_vertex_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # variational averaged

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_vertex_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # variational averaged
    norm  = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):
            norm = max(norm, math.fabs(numerical[iVertex] - analytical[iVertex]))

    return norm

#--------------------------------------------------------

def get_norm_strain_vertex_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # variational averaged

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    kiteAreasOnVertex = fileGrid.variables["kiteAreasOnVertex"][:]
    areaTriangle = fileGrid.variables["areaTriangle"][:]

    verticesOnCell[:] = verticesOnCell[:] - 1
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    fileGrid.close()

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    #strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    #strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    strain11var = fileMPAS.variables["strain11var"][0,:]
    #strain22var = fileMPAS.variables["strain22var"][0,:]
    #strain12var = fileMPAS.variables["strain12var"][0,:]

    fileMPAS.close()

    strain11varAvg = np.zeros(nVertices)
    #strain22varAvg = np.zeros(nVertices)
    #strain12varAvg = np.zeros(nVertices)
    denom = np.zeros(nVertices)
    for iVertex in range(0,nVertices):
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex2 = verticesOnCell[iCell,iVertexOnCell]
                if (iVertex2 == iVertex):
                    strain11varAvg[iVertex] += strain11var[iCell,iVertexOnCell]
                    #strain22varAvg[iVertex] += strain22var[iCell,iVertexOnCell]
                    #strain12varAvg[iVertex] += strain12var[iCell,iVertexOnCell]
                    denom[iVertex] += 1.0
    for iVertex in range(0,nVertices):
        if (denom[iVertex] > 0.0):
            strain11varAvg[iVertex] /= denom[iVertex]
            #strain22varAvg[iVertex] /= denom[iVertex]
            #strain12varAvg[iVertex] /= denom[iVertex]

    if (normType == "l2"):
        normE11 = L2_norm_strain_vertex_avg(strain11varAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_vertex_avg(strain22varAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_vertex_avg(strain12varAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_vertex_avg(strain11varAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_vertex_avg(strain22varAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_vertex_avg(strain12varAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_strain_cell_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # weak averaged

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_cell_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # weak averaged
    norm  = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):
            norm = max(norm, math.fabs(numerical[iVertex] - analytical[iVertex]))

    return norm

#--------------------------------------------------------

def get_norm_strain_cell_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # weak averaged

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]
    kiteAreasOnVertex = fileGrid.variables["kiteAreasOnVertex"][:]
    areaCell = fileGrid.variables["areaCell"][:]
    areaTriangle = fileGrid.variables["areaTriangle"][:]

    verticesOnCell[:] = verticesOnCell[:] - 1
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    fileGrid.close()

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    #strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    #strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    #strain22weak = fileMPAS.variables["strain22weak"][0,:]
    #strain12weak = fileMPAS.variables["strain12weak"][0,:]

    fileMPAS.close()

    strain11weakAvg = np.zeros(nVertices)
    #strain22weakAvg = np.zeros(nVertices)
    #strain12weakAvg = np.zeros(nVertices)
    denom = np.zeros(nVertices)
    for iVertex in range(0,nVertices):
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            if (iCell > 0 and iCell < nCells):
                strain11weakAvg[iVertex] += strain11weak[iCell]
                #strain22weakAvg[iVertex] += strain22weak[iCell]
                #strain12weakAvg[iVertex] += strain12weak[iCell]
                denom[iVertex] += 1.0
    for iVertex in range(0,nVertices):
        if (denom[iVertex] > 0.0):
            strain11weakAvg[iVertex] /= denom[iVertex]
            #strain22weakAvg[iVertex] /= denom[iVertex]
            #strain12weakAvg[iVertex] /= denom[iVertex]

    if (normType == "l2"):
        normE11 = L2_norm_strain_cell_avg(strain11weakAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_cell_avg(strain22weakAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_cell_avg(strain12weakAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_cell_avg(strain11weakAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_cell_avg(strain22weakAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_cell_avg(strain12weakAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    return normE11, normE22, normE12

#--------------------------------------------------------
# stress divergence
#--------------------------------------------------------

def L2_norm_stress_divergence(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_stress_divergence(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):
            norm = max(norm, math.fabs(numerical[iVertex] - analytical[iVertex]))

    return norm

#--------------------------------------------------------

def get_norm_stress_divergence(filenameIC, filename, latitudeLimit, normType):

    fileIC = Dataset(filenameIC, "r")

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    #stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    latVertex = fileMPAS.variables["latVertex"][:]

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    stressDivergenceU = fileMPAS.variables["stressDivergenceU"][0,:]
    #stressDivergenceV = fileMPAS.variables["stressDivergenceV"][0,:]

    if (normType == "l2"):
        normU = L2_norm_stress_divergence(stressDivergenceU, stressDivergenceUAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normV = 0.0#L2_norm_stress_divergence(stressDivergenceV, stressDivergenceVAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    elif (normType == "linf"):
        normU = Linf_norm_stress_divergence(stressDivergenceU, stressDivergenceUAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normV = 0.0#Linf_norm_stress_divergence(stressDivergenceV, stressDivergenceVAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

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

def strain_scaling(axes, normType, label, xlabel, xMin, xMax, yMin):

    resolutions = [2562,10242,40962,163842]

    methods = ["wachspress","pwl","weak","wachspress_avg","pwl_avg","weak_avg"]

    latitudeLimit = 20.0

    dirname = {"wachspress":"wachspress",
               "pwl":"pwl",
               "weak":"weak",
               "wachspress_avg":"wachspress",
               "pwl_avg":"pwl",
               "weak_avg":"weak"}

    # plot options
    legendLabels = {"wachspress":"Wachs.",
                    "pwl":"PWL",
                    "weak":"Weak",
                    "wachspress_avg":"Wachs. Avg.",
                    "pwl_avg":"PWL Avg.",
                    "weak_avg":"Weak Avg."}

    lineColors = {"wachspress":"black",
                  "pwl":"grey",
                  "weak":"red",
                  "wachspress_avg":"blue",
                  "pwl_avg":"green",
                  "weak_avg":"darkturquoise"}

    lineStyles = {"wachspress":"solid",
                  "pwl":"solid",
                  "weak":"solid",
                  "wachspress_avg":"solid",
                  "pwl_avg":"solid",
                  "weak_avg":"solid"}

    markers = {"wachspress":"+",
               "pwl":"x",
               "weak":"^",
               "wachspress_avg":"+",
               "pwl_avg":"x",
               "weak_avg":"^"}

    # plot
    # scaling lines
    scaling_lines(axes, xMin, xMax, yMin)

    plotHandles = []

    for method in methods:

        x = []
        y = []

        for resolution in resolutions:

            filenameGrid = "./strain/grid.%i.nc" %(resolution)
            filenameIC = "./strain/ic_%i.nc" %(resolution)
            filename = "./strain/output_%s_%i/output.2000.nc" %(dirname[method],resolution)

            print(filename, filenameIC)

            if (method == "wachspress" or
                method == "pwl"):
                normE11, normE22, normE12 = get_norm_strain_vertex(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "weak"):
                normE11, normE22, normE12 = get_norm_strain_cell(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "wachspress_avg" or
                  method == "pwl_avg"):
                normE11, normE22, normE12 = get_norm_strain_vertex_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "weak_avg"):
                normE11, normE22, normE12 = get_norm_strain_cell_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)

            x.append(get_resolution(filename, latitudeLimit))
            y.append(normE11)

        plotHandle, = axes.loglog(x,y, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])
        plotHandles.append(plotHandle)

    #axes.legend(frameon=False, loc=2, fontsize=8, handlelength=4)
    legend1 = axes.legend(handles=plotHandles[0:3], frameon=False, loc=2, fontsize=8, handlelength=2)
    legend2 = axes.legend(handles=plotHandles[3:6], frameon=False, loc=4, fontsize=8, handlelength=2)
    axes.add_artist(legend1)

    axes.set_xlabel(xlabel)
    if (normType == "l2"):
        axes.set_ylabel(r"$L_2$ error norm")
    elif (normType == "linf"):
        axes.set_ylabel(r"$L_\infty$ error norm")
    #axes.set_xlim([xMin, xMax])
    ###axes.set_ylim([None, 2e-1])
    axes.set_title(r'%s $\dot{\epsilon}_{11}$' %(label) ,loc="left")
    axes.set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes.set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes.set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes.set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)
    #axes.set_xticklabels(labels=[1e-2,5e-2],minor=False)

#--------------------------------------------------------

def stress_divergence_scaling(axes, normType, label, xlabel, xMin, xMax, yMin):

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
    # scaling lines
    scaling_lines(axes, xMin, xMax, yMin)

    for method in methods:

        x = []
        y = []

        for resolution in resolutions:

            filename = "./stress_divergence/output_%s_%i/output.2000.nc" %(method,resolution)
            filenameIC = "./stress_divergence/ic_%i.nc" %(resolution)

            print(filename, filenameIC)

            normU, normV = get_norm_stress_divergence(filenameIC, filename, latitudeLimit, normType)

            x.append(get_resolution(filename, latitudeLimit))
            y.append(normU)

        axes.loglog(x,y, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])

    axes.legend(frameon=False, loc=4, fontsize=8, handlelength=4)

    axes.set_xlabel(xlabel)
    #axes.set_ylabel(r"$L_2$ error norm")
    axes.set_ylabel(None)
    #axes.set_xlim([8e-3, 8.5e-2])
    axes.set_title(r'%s $(\nabla \cdot \sigma)_u$' %(label), loc="left")
    axes.set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes.set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes.set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes.set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)
    #axes.set_xticklabels(labels=[1e-2,5e-2],minor=False)

#-------------------------------------------------------------------------------

def spherical_operators_scaling():

    cm = 1/2.54  # centimeters in inches
    plt.rc('font',family="Times New Roman")
    plt.rc('mathtext',fontset="stix")
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig, axes = plt.subplots(1, 2, figsize=(15*cm,7*cm))

    strain_scaling(axes[0], "l2", "(a)", "Grid resolution")

    stress_divergence_scaling(axes[1], "l2", "(b)", "Grid resolution")

    #plt.tight_layout()
    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    plt.savefig("spherical_operators_scaling.png",dpi=300)
    plt.savefig("spherical_operators_scaling.eps")

#-------------------------------------------------------------------------------

def spherical_operators_scaling_L2_Linf():

    cm = 1/2.54  # centimeters in inches
    plt.rc('font',family="Times New Roman")
    plt.rc('mathtext',fontset="stix")
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig, axes = plt.subplots(2, 2, figsize=(15*cm,14*cm))

    xMin = 1.5e-2
    xMax = 3e-2
    yMin = 1.5e-4

    strain_scaling(axes[0,0], "l2", "(a)", None, xMin, xMax, yMin)

    xMin = 1e-2
    xMax = 2e-2
    yMin = 7e-3

    stress_divergence_scaling(axes[0,1], "l2", "(b)", None, xMin, xMax, yMin)

    xMin = 1.5e-2
    xMax = 3e-2
    yMin = 1e-3

    strain_scaling(axes[1,0], "linf", "(c)", "Grid resolution", xMin, xMax, yMin)

    xMin = 1e-2
    xMax = 2e-2
    yMin = 2.2e-1

    stress_divergence_scaling(axes[1,1], "linf", "(d)", "Grid resolution", xMin, xMax, yMin)

    #plt.tight_layout()
    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    plt.savefig("spherical_operators_scaling_l2linf.png",dpi=300)
    plt.savefig("spherical_operators_scaling_l2linf.eps")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    spherical_operators_scaling_L2_Linf()
