from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

#---------------------------------------------------------------

def ortho_projection(x, y, z):

    xp = x
    yp = y
    zp = z

    #xp = y
    #yp = z
    #zp = x

    return xp, yp, zp

#---------------------------------------------------------------

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    minval =  1.0e30
    maxval = -1.0e30

    patches = []
    colours = []
    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > math.radians(20.0)):

            useVertex = False
            polygonVertices = []
            for iCellOnVertex in range(0,vertexDegree[iVertex]):

                iCell = cellsOnVertex[iVertex,iCellOnVertex]
                x, y, z = ortho_projection(xCell[iCell], yCell[iCell], zCell[iCell])

                if (x >= minX and x <= maxX and \
                    y >= minY and y <= maxY and \
                    z > 0.0):
                    useVertex = True

                polygonVertices.append((x,y))

            # create patch and add to collection
            if (useVertex):
                polygon = Polygon(polygonVertices)
                patches.append(polygon)

                colours.append(mpasArray[iVertex])

                minval = min(minval,mpasArray[iVertex])
                maxval = max(maxval,mpasArray[iVertex])

    patchCollection = PatchCollection(patches, cmap=cmap, rasterized=True)
    patchCollection.set_array(np.array(colours))
    patchCollection.set_linewidth(0)

    patchCollection.set_clim(vmin=vmin,vmax=vmax)

    return patchCollection, minval, maxval

#---------------------------------------------------------------

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False, cbticks=None, cbticklabels=None):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, array, colourMap, vmin, vmax, minX, maxX, minY, maxY)
    axes.add_collection(patchCollection)
    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

    if (title != None):
        axes.set_title(title, fontsize=8, loc="left")

    if (subfigureLabel != None):
        axes.text(0.02, 0.89, subfigureLabel, verticalalignment='bottom', horizontalalignment='left',transform=axes.transAxes, fontsize=8)

    if (colorbar):
        divider = make_axes_locatable(axes)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        cax = divider.new_horizontal(size="7%",pad=0.05)
        fig.add_axes(cax)
        cb = fig.colorbar(patchCollection,cax=cax)
        if (cbticks is not None):
            cb.set_ticks(cbticks)
            cb.set_ticklabels(cbticklabels)
        if (unityBar):
            cb.ax.set_yticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
        if (sciNote):
            cb.formatter.set_powerlimits((0, 0))
            cb.update_ticks()

#---------------------------------------------------------------

def spherical_operators_map():

    # grid
    fileGrid = Dataset("./strain/x1.40962.grid.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    edgesOnCell = fileGrid.variables["edgesOnCell"][:]

    cellsOnVertex[:] = cellsOnVertex[:] - 1
    verticesOnCell[:] = verticesOnCell[:] - 1
    edgesOnCell[:] = edgesOnCell[:] - 1

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]
    latEdge = fileGrid.variables["latEdge"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    xEdge = fileGrid.variables["xEdge"][:]
    yEdge = fileGrid.variables["yEdge"][:]
    zEdge = fileGrid.variables["zEdge"][:]

    fileGrid.close()

    # variational strain cells
    xVar = []
    yVar = []
    zVar = []
    latVar = []

    nEdgesOnCellVar = []
    verticesOnCellVar = []

    iVar = 0
    nCellsVar = 0
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

            iVertex = verticesOnCell[iCell,iVertexOnCell]
            iEdgeOnCell1 = iVertexOnCell
            iEdgeOnCell2 = iVertexOnCell - 1
            if (iEdgeOnCell2 < 0): iEdgeOnCell2 = nEdgesOnCell[iCell]-1
            iEdge1 = edgesOnCell[iCell,iEdgeOnCell1]
            iEdge2 = edgesOnCell[iCell,iEdgeOnCell2]

            verticesOnCellVar.append([iVar, iVar+1, iVar+2, iVar+3])

            xVar.append(xCell[iCell])
            yVar.append(yCell[iCell])
            zVar.append(zCell[iCell])
            iVar += 1

            xVar.append(xEdge[iEdge1])
            yVar.append(yEdge[iEdge1])
            zVar.append(zEdge[iEdge1])
            iVar += 1

            xVar.append(xVertex[iVertex])
            yVar.append(yVertex[iVertex])
            zVar.append(zVertex[iVertex])
            iVar += 1

            xVar.append(xEdge[iEdge2])
            yVar.append(yEdge[iEdge2])
            zVar.append(zEdge[iEdge2])
            iVar += 1

            latVar.append(0.25 * (latCell[iCell] + latEdge[iEdge1] + latVertex[iVertex] + latEdge[iEdge2]))

            nEdgesOnCellVar.append(4)

            nCellsVar += 1

    xVar = np.array(xVar)
    yVar = np.array(yVar)
    zVar = np.array(zVar)
    latVar = np.array(latVar)
    nEdgesOnCellVar = np.array(nEdgesOnCellVar)
    verticesOnCellVar = np.array(verticesOnCellVar)

    #-----------------------------------------------------------------
    # strains
    #-----------------------------------------------------------------
    print("Strains")

    # ic
    fileIC = Dataset("./strain/ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    strain11VarAnalytical = np.zeros((nCells, maxEdges))
    strain22VarAnalytical = np.zeros((nCells, maxEdges))
    strain12VarAnalytical = np.zeros((nCells, maxEdges))
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            strain11VarAnalytical[iCell,iVertexOnCell] = strain11VertexAnalytical[iVertex]
            strain22VarAnalytical[iCell,iVertexOnCell] = strain22VertexAnalytical[iVertex]
            strain12VarAnalytical[iCell,iVertexOnCell] = strain12VertexAnalytical[iVertex]

    print("  Analytical: ",
          np.amin(strain11VertexAnalytical), np.amax(strain11VertexAnalytical),
          np.amin(strain22VertexAnalytical), np.amax(strain22VertexAnalytical),
          np.amin(strain12VertexAnalytical), np.amax(strain12VertexAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./strain/output_wachspress_40962/output.2000.nc","r")

    strain11varWachspress = fileWach.variables["strain11var"][0,:]
    strain22varWachspress = fileWach.variables["strain22var"][0,:]
    strain12varWachspress = fileWach.variables["strain12var"][0,:]

    strain11varWachspressDiff = []
    strain22varWachspressDiff = []
    strain12varWachspressDiff = []
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            strain11varWachspressDiff.append(strain11varWachspress[iCell,iVertexOnCell] - strain11VarAnalytical[iCell,iVertexOnCell])
            strain22varWachspressDiff.append(strain22varWachspress[iCell,iVertexOnCell] - strain22VarAnalytical[iCell,iVertexOnCell])
            strain12varWachspressDiff.append(strain12varWachspress[iCell,iVertexOnCell] - strain12VarAnalytical[iCell,iVertexOnCell])

    print("  Wachs: ",
          np.amin(strain11varWachspressDiff), np.amax(strain11varWachspressDiff),
          np.amin(strain22varWachspressDiff), np.amax(strain22varWachspressDiff),
          np.amin(strain12varWachspressDiff), np.amax(strain12varWachspressDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./strain/output_pwl_40962/output.2000.nc","r")

    strain11varPWL = filePWL.variables["strain11var"][0,:]
    strain22varPWL = filePWL.variables["strain22var"][0,:]
    strain12varPWL = filePWL.variables["strain12var"][0,:]

    strain11varPWLDiff = []
    strain22varPWLDiff = []
    strain12varPWLDiff = []
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            strain11varPWLDiff.append(strain11varPWL[iCell,iVertexOnCell] - strain11VarAnalytical[iCell,iVertexOnCell])
            strain22varPWLDiff.append(strain22varPWL[iCell,iVertexOnCell] - strain22VarAnalytical[iCell,iVertexOnCell])
            strain12varPWLDiff.append(strain12varPWL[iCell,iVertexOnCell] - strain12VarAnalytical[iCell,iVertexOnCell])

    print("  PWL:   ",
          np.amin(strain11varPWLDiff), np.amax(strain11varPWLDiff),
          np.amin(strain22varPWLDiff), np.amax(strain22varPWLDiff),
          np.amin(strain12varPWLDiff), np.amax(strain12varPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./strain/output_weak_40962/output.2000.nc","r")

    strain11weakWeak = fileWeak.variables["strain11weak"][0,:]
    strain22weakWeak = fileWeak.variables["strain22weak"][0,:]
    strain12weakWeak = fileWeak.variables["strain12weak"][0,:]

    strain11weakWeakDiff = strain11weakWeak - strain11CellAnalytical
    strain22weakWeakDiff = strain22weakWeak - strain22CellAnalytical
    strain12weakWeakDiff = strain12weakWeak - strain12CellAnalytical

    print("  Weak:  ",
          np.amin(strain11weakWeakDiff), np.amax(strain11weakWeakDiff),
          np.amin(strain22weakWeakDiff), np.amax(strain22weakWeakDiff),
          np.amin(strain12weakWeakDiff), np.amax(strain12weakWeakDiff))

    fileWeak.close()

    # Weak Wachs
    fileWeakWachs = Dataset("./strain/output_weakwachs_40962/output.2000.nc","r")

    strain11weakWachs = fileWeakWachs.variables["strain11varAvg"][0,:]
    strain22weakWachs = fileWeakWachs.variables["strain22varAvg"][0,:]
    strain12weakWachs = fileWeakWachs.variables["strain12varAvg"][0,:]

    strain11weakWachsDiff = strain11weakWachs - strain11VertexAnalytical
    strain22weakWachsDiff = strain22weakWachs - strain22VertexAnalytical
    strain12weakWachsDiff = strain12weakWachs - strain12VertexAnalytical

    print("  WeakWachs:  ",
          np.amin(strain11weakWachsDiff), np.amax(strain11weakWachsDiff),
          np.amin(strain22weakWachsDiff), np.amax(strain22weakWachsDiff),
          np.amin(strain12weakWachsDiff), np.amax(strain12weakWachsDiff))

    fileWeakWachs.close()

    #-----------------------------------------------------------------
    # divergence stress
    #-----------------------------------------------------------------
    print("Stress divergence")

    # ic
    fileIC = Dataset("./stress_divergence/ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    print("  Stress divergence: ",
          np.amin(stressDivergenceUAnalytical), np.amax(stressDivergenceUAnalytical),
          np.amin(stressDivergenceVAnalytical), np.amax(stressDivergenceVAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./stress_divergence/output_wachspress_40962/output.2000.nc","r")

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("  Wachs: ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    # Wachspress alt
    fileWach = Dataset("./stress_divergence/output_wachspress_alt_40962/output.2000.nc","r")

    stressDivergenceUWachAlt = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWachAlt = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachAltDiff = (stressDivergenceUWachAlt - stressDivergenceUAnalytical)
    stressDivergenceVWachAltDiff = (stressDivergenceVWachAlt - stressDivergenceVAnalytical)

    print("  Wachs: ",
          np.amin(stressDivergenceUWachAltDiff), np.amax(stressDivergenceUWachAltDiff),
          np.amin(stressDivergenceVWachAltDiff), np.amax(stressDivergenceVWachAltDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./stress_divergence/output_pwl_40962/output.2000.nc","r")

    stressDivergenceUPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLDiff = (stressDivergenceUPWL - stressDivergenceUAnalytical)
    stressDivergenceVPWLDiff = (stressDivergenceVPWL - stressDivergenceVAnalytical)

    print("  PWL:   ",
          np.amin(stressDivergenceUPWLDiff), np.amax(stressDivergenceUPWLDiff),
          np.amin(stressDivergenceVPWLDiff), np.amax(stressDivergenceVPWLDiff))

    filePWL.close()

    # PWL Alt
    filePWL = Dataset("./stress_divergence/output_pwl_alt_40962/output.2000.nc","r")

    stressDivergenceUPWLAlt = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWLAlt = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLAltDiff = (stressDivergenceUPWLAlt - stressDivergenceUAnalytical)
    stressDivergenceVPWLAltDiff = (stressDivergenceVPWLAlt - stressDivergenceVAnalytical)

    print("  PWL:   ",
          np.amin(stressDivergenceUPWLAltDiff), np.amax(stressDivergenceUPWLAltDiff),
          np.amin(stressDivergenceVPWLAltDiff), np.amax(stressDivergenceVPWLAltDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./stress_divergence/output_weak_40962/output.2000.nc","r")

    stressDivergenceUWeak = fileWeak.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeak = fileWeak.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakDiff = (stressDivergenceUWeak - stressDivergenceUAnalytical)
    stressDivergenceVWeakDiff = (stressDivergenceVWeak - stressDivergenceVAnalytical)

    print("  Weak:  ",
          np.amin(stressDivergenceUWeakDiff), np.amax(stressDivergenceUWeakDiff),
          np.amin(stressDivergenceVWeakDiff), np.amax(stressDivergenceVWeakDiff))

    fileWeak.close()


    cm = 1/2.54  # centimeters in inches
    plt.rc('font', family="Times New Roman")
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


    fig, axes = plt.subplots(5,4,figsize=(15*cm,17.5*cm))

    xDetailMin = -0.1
    xDetailMax =  0.1
    yDetailMin = -0.1
    yDetailMax =  0.1

    minStrain = -3.4
    maxStrain =  3.4

    minStrainDiff = -0.12
    maxStrainDiff =  0.12

    minStrainWeakDiff = -0.004
    maxStrainWeakDiff =  0.004

    minStressDiv = -20.0
    maxStressDiv =  20.0

    minStressDivDiff = -0.8
    maxStressDivDiff =  0.8

    # Velocities/
    firstRowColorbars = True
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(a) $u^\prime$', None, False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(b) $v^\prime$', None, firstRowColorbars)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(c) $(\nabla \cdot \sigma)_{u^\prime}$ Analytical', None, False)
    plot_subfigure(axes[0,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(d) $(\nabla \cdot \sigma)_{v^\prime}$ Analytical', None, firstRowColorbars)

    # Analytical strains
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(e) $\dot{\epsilon}_{11}$', None, False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(f) $\dot{\epsilon}_{22}$', None, False)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(g) $\dot{\epsilon}_{12}$', None, True)
    axes[1,3].axis('off')

    # Wachspress
    plot_subfigure(axes[2,0], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(h) $\dot{\epsilon}_{11}$ Wachs.', None, False)
    plot_subfigure(axes[2,1], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(i) $\dot{\epsilon}_{11}$ Wachs.', None, True)
    #plot_subfigure(axes[2,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
    #               False, False, r'(h) $\dot{\epsilon}_{11}$ Wachs.', None, False)
    #plot_subfigure(axes[2,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
    #               False, False, r'(i) $\dot{epsilon}_{11}$ Wachs.', None, True)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(j) $(\nabla \cdot \sigma)_{u^\prime}$ Wachs.', None, False)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(k) $(\nabla \cdot \sigma)_{u^\prime}$ Wachs.', None, True, cbticks=[-0.8,-0.4,0.0,0.4,0.8], cbticklabels=["-0.8","-0.4","0.0","0.4","0.8"])

    # PWL
    plot_subfigure(axes[3,0], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varPWLDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(l) $\dot{\epsilon}_{11}$ PWL', None, False)
    plot_subfigure(axes[3,1], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varPWLDiff, minStrainDiff, maxStrainDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(m) $\dot{\epsilon}_{11}$ PWL', None, True)
    #plot_subfigure(axes[3,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varPWLDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
    #               False, False, r'(l) $\dot{\epsilon}_{11}$ PWL', None, False)
    #plot_subfigure(axes[3,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varPWLDiff, minStrainDiff, maxStrainDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
    #               False, False, r'(m) $\dot{\epsilon}_{11}$ PWL', None, True)
    plot_subfigure(axes[3,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(n) $(\nabla \cdot \sigma)_{u^\prime}$ PWL', None, False)
    plot_subfigure(axes[3,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLDiff, minStressDivDiff, maxStressDivDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(o) $(\nabla \cdot \sigma)_{u^\prime}$ PWL', None, True, cbticks=[-0.8,-0.4,0.0,0.4,0.8], cbticklabels=["-0.8","-0.4","0.0","0.4","0.8"])

    # Weak
    plot_subfigure(axes[4,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain11weakWeakDiff[:]*1e2, minStrainWeakDiff*1e2, maxStrainWeakDiff*1e2, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(p) $\dot{\epsilon}_{11}$ Weak ($\times 10^{-2}$)', None, False)
    plot_subfigure(axes[4,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain11weakWeakDiff[:]*1e2, minStrainWeakDiff*1e2, maxStrainWeakDiff*1e2, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(q) $\dot{\epsilon}_{11}$ Weak ($\times 10^{-2}$)', None, True)
    plot_subfigure(axes[4,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(r) $(\nabla \cdot \sigma)_{u^\prime}$ Weak', None, False)
    plot_subfigure(axes[4,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(s) $(\nabla \cdot \sigma)_{u^\prime}$ Weak', None, True, cbticks=[-0.8,-0.4,0.0,0.4,0.8], cbticklabels=["-0.8","-0.4","0.0","0.4","0.8"])



    #plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0)

    plt.subplots_adjust(left=0.015,
                        bottom=0.015,
                        right=0.925,
                        top=0.97,
                        wspace=0.43,
                        hspace=0.2)

    positionTypeWrite = "read"

    if (positionTypeWrite == "read"):

        # vertical
        filein = open("axis_position_vertical.txt","r")
        lines = filein.readlines()
        filein.close()
        y0 = []
        y1 = []
        for i in range(0,5):
            y0.append(float(lines[i].split()[1]))
            y1.append(float(lines[i].split()[2]))

        # horizontal
        filein = open("axis_position_horizontal.txt","r")
        lines = filein.readlines()
        filein.close()
        x0 = []
        x1 = []
        for j in range(0,4):
            x0.append(float(lines[j].split()[1]))
            x1.append(float(lines[j].split()[2]))

        cbs = [(0,1),(2,1),(3,1),(4,1),(1,2),(0,3),(2,3),(3,3),(4,3)]
        for i in range(0,5):
            for j in range(0,4):
                if ((i,j) in cbs):
                    width = (x1[j]-x0[j]) * 1.12
                else:
                    width = x1[j]-x0[j]
                pos = axes[i,j].get_position(original=False)
                newPos = [x0[j], y0[i], width, y1[i]-y0[i]]
                axes[i,j].set_position(newPos) # set a new position



    #plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0)

    if (positionTypeWrite == "write"):

        # vertical
        fileout = open("axis_position_vertical.txt","w")
        for i in range(0,5):
            pos = axes[i,0].get_position(original=False)
            fileout.write("%i %f %f\n" %(i,pos.y0,pos.y1))
        fileout.close()

        # horizontal
        fileout = open("axis_position_horizontal.txt","w")
        for j in range(0,4):
            pos = axes[0,j].get_position(original=False)
            fileout.write("%i %f %f\n" %(i,pos.x0,pos.x1))
        fileout.close()



    for i in range(0,5):
        for j in range(0,4):
            pos = axes[i,j].get_position(original=False)
            print(i,j,pos)




    plt.savefig("spherical_operators_map.png",dpi=300)
    plt.savefig("spherical_operators_map.eps")

    plt.clf()
    plt.cla()
    plt.close()



    # alt vs orig
    fig, axes = plt.subplots(2,2,figsize=(8*cm,7.5*cm))

    minStressDivDiff = -0.6
    maxStressDivDiff =  0.6

    # orig
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachAltDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(a) $(\nabla \cdot \sigma)_{u^\prime}$ Wachs. Alt.', None, False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachAltDiff, minStressDivDiff, maxStressDivDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(b) $(\nabla \cdot \sigma)_{u^\prime}$ Wachs. Alt.', None, True, cbticks=[-0.6,-0.3,0.0,0.3,0.6], cbticklabels=["-0.6","-0.3","0.0","0.3","0.6"])

    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLAltDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'(c) $(\nabla \cdot \sigma)_{u^\prime}$ PWL Alt.', None, False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLAltDiff, minStressDivDiff, maxStressDivDiff, xDetailMin, xDetailMax, yDetailMin, yDetailMax, \
                   False, False, r'(d) $(\nabla \cdot \sigma)_{u^\prime}$ PWL Alt.', None, True, cbticks=[-0.6,-0.3,0.0,0.3,0.6], cbticklabels=["-0.6","-0.3","0.0","0.3","0.6"])

    plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0)
    plt.savefig("spherical_operators_map_comp.png",dpi=300)
    plt.savefig("spherical_operators_map_comp.eps")

    plt.clf()
    plt.cla()
    plt.close()


#---------------------------------------------------------------

if __name__ == "__main__":

    spherical_operators_map()
