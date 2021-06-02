from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#---------------------------------------------------------------

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        if (latVertex[iVertex] > math.radians(20.0)):

            polygonVertices = []

            useVertex = False
            for iCellOnVertex in range(0,vertexDegree[iVertex]):

                iCell = cellsOnVertex[iVertex,iCellOnVertex]

                polygonVertices.append((xCell[iCell],yCell[iCell]))

                if (xCell[iCell] >= minX and xCell[iCell] <= maxX and \
                    yCell[iCell] >= minY and yCell[iCell] <= maxY):
                    useVertex = True

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

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

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
        axes.set_title(title, fontsize=8)

    if (subfigureLabel != None):
        axes.text(0.02, 0.89, subfigureLabel, verticalalignment='bottom', horizontalalignment='left',transform=axes.transAxes, fontsize=8)

    if (colorbar):
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(patchCollection,cax=cax)
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
    fileWach = Dataset("./stress_divergence/output_wachspress_alt_40962/output.2000.nc","r")

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("  Wachs: ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./stress_divergence/output_pwl_alt_40962/output.2000.nc","r")

    stressDivergenceUPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLDiff = (stressDivergenceUPWL - stressDivergenceUAnalytical)
    stressDivergenceVPWLDiff = (stressDivergenceVPWL - stressDivergenceVAnalytical)

    print("  PWL:   ",
          np.amin(stressDivergenceUPWLDiff), np.amax(stressDivergenceUPWLDiff),
          np.amin(stressDivergenceVPWLDiff), np.amax(stressDivergenceVPWLDiff))

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




    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5


    fig, axes = plt.subplots(5, 4)
    fig.set_size_inches(7.5, 8)

    minStrain = -3.4
    maxStrain =  3.4

    minStrainDiff = -0.12
    maxStrainDiff =  0.12

    minStrainWeakDiff = -0.005
    maxStrainWeakDiff =  0.005

    minStressDiv = -20.0
    maxStressDiv =  20.0

    minStressDivDiff = -0.2
    maxStressDivDiff =  0.2

    # Velocities/
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(b)', True)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Analytical', '(c)', False)
    plot_subfigure(axes[0,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{v^\prime}$ Analytical', '(d)', True)

    # Analytical strains
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$', '(e)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$', '(f)', False)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$', '(g)', True)
    axes[1,3].axis('off')

    # Wachspress
    plot_subfigure(axes[2,0], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ Wachs.', '(h)', False)
    plot_subfigure(axes[2,1], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ Wachs.', '(i)', True)
    #plot_subfigure(axes[2,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
    #               False, False, r'$\epsilon_{11}$ Wachs.', '(h)', False)
    #plot_subfigure(axes[2,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varWachspressDiff, minStrainDiff, maxStrainDiff, -0.2, 0.2, -0.2, 0.2, \
    #               False, False, r'$\epsilon_{11}$ Wachs.', '(i)', True)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Wachs.', '(j)', False)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Wachs.', '(k)', True)

    # PWL
    plot_subfigure(axes[3,0], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varPWLDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ PWL', '(l)', False)
    plot_subfigure(axes[3,1], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, latVar, strain11varPWLDiff, minStrainDiff, maxStrainDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ PWL', '(m)', True)
    #plot_subfigure(axes[3,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varPWLDiff, minStrainDiff, maxStrainDiff, -1.0, 1.0, -1.0, 1.0, \
    #               False, False, r'$\epsilon_{11}$ PWL', '(l)', False)
    #plot_subfigure(axes[3,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varPWLDiff, minStrainDiff, maxStrainDiff, -0.2, 0.2, -0.2, 0.2, \
    #               False, False, r'$\epsilon_{11}$ PWL', '(m)', True)
    plot_subfigure(axes[3,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ PWL', '(n)', False)
    plot_subfigure(axes[3,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ PWL', '(o)', True)

    # Weak
    plot_subfigure(axes[4,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain11weakWeakDiff, minStrainWeakDiff, maxStrainWeakDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ Weak', '(p)', False)
    plot_subfigure(axes[4,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain11weakWeakDiff, minStrainWeakDiff, maxStrainWeakDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ Weak', '(q)', True)
    plot_subfigure(axes[4,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Weak', '(r)', False)
    plot_subfigure(axes[4,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Weak', '(s)', True)



    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("spherical_operators_map.png",dpi=400)

    plt.clf()
    plt.cla()
    plt.close()

#---------------------------------------------------------------

if __name__ == "__main__":

    spherical_operators_map()
