from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

degreesToRadians = math.pi / 180.0

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#---------------------------------------------------------------

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertexx, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        polygonVertices = []

        useVertex = False
        for iCellOnVertex in range(0,vertexDegree[iVertex]):

            iCell = cellsOnVertex[iVertex,iCellOnVertex]

            polygonVertices.append((xCell[iCell],yCell[iCell]))

            if (xCell[iCell] >= minX and xCell[iCell] <= maxX and \
                yCell[iCell] >= minY and yCell[iCell] <= maxY):
                useVertex = True

        if (useVertex and (useVertexx[iVertex] == 1)):
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

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, array, colourMap, vmin, vmax, minX, maxX, minY, maxY)
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

def strain_map():

    # grid quad
    fileGrid = Dataset("grid_hex_0082x0094.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    cellsOnCell = fileGrid.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    fileGrid.close()

    # ic hex
    fileIC = Dataset("ic_hex_0082x0094.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    print("Stress divergence: ",
          np.amin(strain11VertexAnalytical), np.amax(strain11VertexAnalytical),
          np.amin(strain22VertexAnalytical), np.amax(strain22VertexAnalytical),
          np.amin(strain12VertexAnalytical), np.amax(strain12VertexAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_hex_wachspress_0082x0094/output.2000.nc","r")

    interiorCell = fileWach.variables["interiorCell"][0,:]

    strain11varAvgWachs = fileWach.variables["strain11varAvg"][0,:]
    strain22varAvgWachs = fileWach.variables["strain22varAvg"][0,:]
    strain12varAvgWachs = fileWach.variables["strain12varAvg"][0,:]

    strain11varAvgWachsDiff = (strain11varAvgWachs - strain11VertexAnalytical)
    strain22varAvgWachsDiff = (strain22varAvgWachs - strain22VertexAnalytical)
    strain12varAvgWachsDiff = (strain12varAvgWachs - strain12VertexAnalytical)

    print("Wachs: ",
          np.amin(strain11varAvgWachsDiff), np.amax(strain11varAvgWachsDiff),
          np.amin(strain22varAvgWachsDiff), np.amax(strain22varAvgWachsDiff),
          np.amin(strain12varAvgWachsDiff), np.amax(strain12varAvgWachsDiff))

    fileWach.close()

    useVertex = np.ones(nVertices,dtype="i")
    useCell = np.ones(nCells,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            useCell[iCell] = 0
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                useCell[iCell2] = 0
                if (iCell2 < nCells):
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell2]):
                        iVertex = verticesOnCell[iCell2,iVertexOnCell]
                        useVertex[iVertex] = 0

    # PWL
    filePWL = Dataset("./output_hex_pwl_0082x0094/output.2000.nc","r")

    strain11varAvgPWL = filePWL.variables["strain11varAvg"][0,:]
    strain22varAvgPWL = filePWL.variables["strain22varAvg"][0,:]
    strain12varAvgPWL = filePWL.variables["strain12varAvg"][0,:]

    strain11varAvgPWLDiff = (strain11varAvgPWL - strain11VertexAnalytical)
    strain22varAvgPWLDiff = (strain22varAvgPWL - strain22VertexAnalytical)
    strain12varAvgPWLDiff = (strain12varAvgPWL - strain12VertexAnalytical)

    print("PWL: ",
          np.amin(strain11varAvgPWLDiff), np.amax(strain11varAvgPWLDiff),
          np.amin(strain22varAvgPWLDiff), np.amax(strain22varAvgPWLDiff),
          np.amin(strain12varAvgPWLDiff), np.amax(strain12varAvgPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./output_hex_weak_0082x0094/output.2000.nc","r")

    strain11weakWeak = fileWeak.variables["strain11weak"][0,:]
    strain22weakWeak = fileWeak.variables["strain22weak"][0,:]
    strain12weakWeak = fileWeak.variables["strain12weak"][0,:]

    strain11weakWeakDiff = (strain11weakWeak - strain11CellAnalytical)
    strain22weakWeakDiff = (strain22weakWeak - strain22CellAnalytical)
    strain12weakWeakDiff = (strain12weakWeak - strain12CellAnalytical)

    print("Weak: ",
          np.amin(strain11weakWeakDiff), np.amax(strain11weakWeakDiff),
          np.amin(strain22weakWeakDiff), np.amax(strain22weakWeakDiff),
          np.amin(strain12weakWeakDiff), np.amax(strain12weakWeakDiff))

    fileWeak.close()


    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minVelocity = -1.0
    maxVelocity =  1.0

    # sinusoid
    #minStrain = -20.0
    #maxStrain =  20.0

    #minStrainDiff = -2.0
    #maxStrainDiff =  2.0

    # linear
    minStrain = None#-1.0
    maxStrain = None# 1.0

    minStrainDiff = None#-1e-2
    maxStrainDiff = None# 1e-2

    fig, axes = plt.subplots(4, 3)

    fig.set_size_inches(7, 6.75)

    # analytical
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain11VertexAnalytical, minStrain, maxStrain, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ Analytical', r'(a)', True)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain22VertexAnalytical, minStrain, maxStrain, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ Analytical', r'(b)', True)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain12VertexAnalytical, minStrain, maxStrain, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ Analytical', r'(c)', True)

    # Wachspress
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain11varAvgWachs, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ Wachs.', r'(d)', True)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain22varAvgWachs, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ Wachs.', r'(e)', True)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain12varAvgWachs, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ Wachs.', r'(f)', True)

    # PWL
    plot_subfigure(axes[2,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain11varAvgPWL, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ PWL', r'(g)', True)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain22varAvgPWL, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ PWL', r'(h)', True)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain12varAvgPWL, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ PWL', r'(i)', True)

    # Weak
    plot_subfigure(axes[3,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, useCell, strain11weakWeak, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ Weak', r'(j)', True)
    plot_subfigure(axes[3,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, useCell, strain22weakWeak, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ Weak', r'(k)', True)
    plot_subfigure(axes[3,2], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, useCell, strain12weakWeak, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ Weak', r'(l)', True)

    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_map.png",dpi=400)

    plt.clf()
    plt.cla()
    plt.close()



#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_map()
