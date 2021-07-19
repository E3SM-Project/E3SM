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

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        if (latVertex[iVertex] > 20.0*degreesToRadians):

            polygonVertices = []

            useVertex = False
            for iCellOnVertex in range(0,vertexDegree[iVertex]):

                iCell = cellsOnVertex[iVertex,iCellOnVertex] - 1

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

def strain_map():

    # grid
    fileGrid = Dataset("x1.40962.grid.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    # ic
    fileIC = Dataset("ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    print("Strain: ",
          np.amin(strain11VertexAnalytical), np.amax(strain11VertexAnalytical),
          np.amin(strain22VertexAnalytical), np.amax(strain22VertexAnalytical),
          np.amin(strain12VertexAnalytical), np.amax(strain12VertexAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_wachspress_40962/output.2000.nc","r")

    strain11varWachspress = fileWach.variables["strain11varAvg"][0,:]
    strain22varWachspress = fileWach.variables["strain22varAvg"][0,:]
    strain12varWachspress = fileWach.variables["strain12varAvg"][0,:]

    strain11varWachspressDiff = strain11varWachspress - strain11VertexAnalytical
    strain22varWachspressDiff = strain22varWachspress - strain22VertexAnalytical
    strain12varWachspressDiff = strain12varWachspress - strain12VertexAnalytical

    print("Wachs: ",
          np.amin(strain11varWachspressDiff), np.amax(strain11varWachspressDiff),
          np.amin(strain22varWachspressDiff), np.amax(strain22varWachspressDiff),
          np.amin(strain12varWachspressDiff), np.amax(strain12varWachspressDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./output_pwl_40962/output.2000.nc","r")

    strain11varPWL = filePWL.variables["strain11varAvg"][0,:]
    strain22varPWL = filePWL.variables["strain22varAvg"][0,:]
    strain12varPWL = filePWL.variables["strain12varAvg"][0,:]

    strain11varPWLDiff = strain11varPWL - strain11VertexAnalytical
    strain22varPWLDiff = strain22varPWL - strain22VertexAnalytical
    strain12varPWLDiff = strain12varPWL - strain12VertexAnalytical

    print("PWL:   ",
          np.amin(strain11varPWLDiff), np.amax(strain11varPWLDiff),
          np.amin(strain22varPWLDiff), np.amax(strain22varPWLDiff),
          np.amin(strain12varPWLDiff), np.amax(strain12varPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./output_weak_40962/output.2000.nc","r")

    strain11weakWeak = fileWeak.variables["strain11weak"][0,:]
    strain22weakWeak = fileWeak.variables["strain22weak"][0,:]
    strain12weakWeak = fileWeak.variables["strain12weak"][0,:]

    strain11weakWeakDiff = strain11weakWeak - strain11CellAnalytical
    strain22weakWeakDiff = strain22weakWeak - strain22CellAnalytical
    strain12weakWeakDiff = strain12weakWeak - strain12CellAnalytical

    print("Weak:  ",
          np.amin(strain11weakWeakDiff), np.amax(strain11weakWeakDiff),
          np.amin(strain22weakWeakDiff), np.amax(strain22weakWeakDiff),
          np.amin(strain12weakWeakDiff), np.amax(strain12weakWeakDiff))

    fileWeak.close()

    # Weak Wachs
    fileWeakWachs = Dataset("./output_weakwachs_40962/output.2000.nc","r")

    strain11weakWachs = fileWeakWachs.variables["strain11varAvg"][0,:]
    strain22weakWachs = fileWeakWachs.variables["strain22varAvg"][0,:]
    strain12weakWachs = fileWeakWachs.variables["strain12varAvg"][0,:]

    strain11weakWachsDiff = strain11weakWachs - strain11VertexAnalytical
    strain22weakWachsDiff = strain22weakWachs - strain22VertexAnalytical
    strain12weakWachsDiff = strain12weakWachs - strain12VertexAnalytical

    print("WeakWachs:  ",
          np.amin(strain11weakWachsDiff), np.amax(strain11weakWachsDiff),
          np.amin(strain22weakWachsDiff), np.amax(strain22weakWachsDiff),
          np.amin(strain12weakWachsDiff), np.amax(strain12weakWachsDiff))

    fileWeakWachs.close()



    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5


    fig, axes = plt.subplots(6, 6)
    fig.set_size_inches(9, 8)

    minStrain = -3.3
    maxStrain =  3.3

    minDiff = -0.043
    maxDiff =  0.043

    # Velocities
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(b)', True)
    axes[0,2].axis('off')
    axes[0,3].axis('off')
    axes[0,4].axis('off')
    axes[0,5].axis('off')

    # Analytical strains
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$', '(c)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11VertexAnalytical, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$', '(d)', False)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$', '(e)', False)
    plot_subfigure(axes[1,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22VertexAnalytical, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$', '(f)', False)
    plot_subfigure(axes[1,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$', '(g)', False)
    plot_subfigure(axes[1,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12VertexAnalytical, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$', '(h)', True)

    # Wachspress
    plot_subfigure(axes[2,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varWachspressDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ Wachs.', '(i)', False)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varWachspressDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ Wachs.', '(j)', False)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22varWachspressDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$ Wachs.', '(k)', False)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22varWachspressDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$ Wachs.', '(l)', False)
    plot_subfigure(axes[2,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12varWachspressDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$ Wachs.', '(m)', False)
    plot_subfigure(axes[2,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12varWachspressDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$ Wachs.', '(n)', True)

    # PWL
    plot_subfigure(axes[3,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varPWLDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ PWL', '(o)', False)
    plot_subfigure(axes[3,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11varPWLDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ PWL', '(p)', False)
    plot_subfigure(axes[3,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22varPWLDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$ PWL', '(q)', False)
    plot_subfigure(axes[3,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22varPWLDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$ PWL', '(r)', False)
    plot_subfigure(axes[3,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12varPWLDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$ PWL', '(s)', False)
    plot_subfigure(axes[3,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12varPWLDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$ PWL', '(t)', True)

    # Weak
    plot_subfigure(axes[4,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain11weakWeakDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ Weak', '(u)', False)
    plot_subfigure(axes[4,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain11weakWeakDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ Weak', '(v)', False)
    plot_subfigure(axes[4,2], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain22weakWeakDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$ Weak', '(w)', False)
    plot_subfigure(axes[4,3], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain22weakWeakDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$ Weak', '(x)', False)
    plot_subfigure(axes[4,4], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain12weakWeakDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$ Weak', '(y)', False)
    plot_subfigure(axes[4,5], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, latCell, strain12weakWeakDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$ Weak', '(z)', True)

    # Weak Wachs
    plot_subfigure(axes[5,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11weakWachsDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ WeakWachs.', '(aa)', False)
    plot_subfigure(axes[5,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11weakWachsDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ WeakWachs.', '(bb)', False)
    plot_subfigure(axes[5,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22weakWachsDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$ WeakWachs.', '(cc)', False)
    plot_subfigure(axes[5,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22weakWachsDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$ WeakWachs.', '(dd)', False)
    plot_subfigure(axes[5,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12weakWachsDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$ WeakWachs.', '(ee)', False)
    plot_subfigure(axes[5,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12weakWachsDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$ WeakWachs.', '(ff)', True)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_map.png",dpi=400)
    #plt.savefig("strain_map_3.png", bbox_inches="tight",dpi=2000)

    plt.clf()
    plt.cla()
    plt.close()

#---------------------------------------------------------------

if __name__ == "__main__":

    strain_map()
