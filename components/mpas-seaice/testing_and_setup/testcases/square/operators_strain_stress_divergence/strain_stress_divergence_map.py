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

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertexx, latVertex, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        polygonVertices = []

        useVertex = False
        for iCellOnVertex in range(0,vertexDegree):

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

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, array, colourMap, vmin, vmax, minX, maxX, minY, maxY)
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

def strain_stress_divergence_map():

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

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    print("Stress divergence: ",
          np.amin(stressDivergenceUAnalytical), np.amax(stressDivergenceUAnalytical),
          np.amin(stressDivergenceVAnalytical), np.amax(stressDivergenceVAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_hex_wachsavg_0082x0094/output.2000.nc","r")

    interiorCell = fileWach.variables["interiorCell"][0,:]

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("Wachs: ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    useVertex = np.ones(nVertices,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                if (iCell2 < nCells):
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell2]):
                        iVertex = verticesOnCell[iCell2,iVertexOnCell]
                        useVertex[iVertex] = 0

    # PWL
    filePWL = Dataset("./output_hex_pwlavg_0082x0094/output.2000.nc","r")

    stressDivergenceUPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLDiff = (stressDivergenceUPWL - stressDivergenceUAnalytical)
    stressDivergenceVPWLDiff = (stressDivergenceVPWL - stressDivergenceVAnalytical)

    print("PWL:   ",
          np.amin(stressDivergenceUPWLDiff), np.amax(stressDivergenceUPWLDiff),
          np.amin(stressDivergenceVPWLDiff), np.amax(stressDivergenceVPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./output_hex_weak_0082x0094/output.2000.nc","r")

    stressDivergenceUWeak = fileWeak.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeak = fileWeak.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakDiff = (stressDivergenceUWeak - stressDivergenceUAnalytical)
    stressDivergenceVWeakDiff = (stressDivergenceVWeak - stressDivergenceVAnalytical)

    print("Weak:  ",
          np.amin(stressDivergenceUWeakDiff), np.amax(stressDivergenceUWeakDiff),
          np.amin(stressDivergenceVWeakDiff), np.amax(stressDivergenceVWeakDiff))

    fileWeak.close()

    # Weak Wachspress
    fileWach = Dataset("./output_hex_weakwachs_0082x0094/output.2000.nc","r")

    stressDivergenceUWeakWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeakWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakWachDiff = (stressDivergenceUWeakWach - stressDivergenceUAnalytical)
    stressDivergenceVWeakWachDiff = (stressDivergenceVWeakWach - stressDivergenceVAnalytical)

    print("Wachs: ",
          np.amin(stressDivergenceUWeakWachDiff), np.amax(stressDivergenceUWeakWachDiff),
          np.amin(stressDivergenceVWeakWachDiff), np.amax(stressDivergenceVWeakWachDiff))

    fileWach.close()

    # Weak PWL
    filePWL = Dataset("./output_hex_weakpwl_0082x0094/output.2000.nc","r")

    stressDivergenceUWeakPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeakPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakPWLDiff = (stressDivergenceUWeakPWL - stressDivergenceUAnalytical)
    stressDivergenceVWeakPWLDiff = (stressDivergenceVWeakPWL - stressDivergenceVAnalytical)

    print("PWL:   ",
          np.amin(stressDivergenceUWeakPWLDiff), np.amax(stressDivergenceUWeakPWLDiff),
          np.amin(stressDivergenceVWeakPWLDiff), np.amax(stressDivergenceVWeakPWLDiff))

    filePWL.close()


    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minVelocity = -1.0
    maxVelocity =  1.0

    minStressDiv = -750.0
    maxStressDiv =  750.0

    minStressDivDiff = -40.0
    maxStressDivDiff =  40.0


    fig, axes = plt.subplots(4, 4)

    fig.set_size_inches(7, 6.75)

    plot_subfigure(axes[0,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, uVelocity, minVelocity, maxVelocity, xMin, xMax, yMin, yMax, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, vVelocity, minVelocity, maxVelocity, xMin, xMax, yMin, yMax, \
                   False, False, r'$v^\prime$', '(b)', False)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUAnalytical, minStressDiv, maxStressDiv, xMin, xMax, yMin, yMax, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Analytical', r'(c)$\times20$', False)
    plot_subfigure(axes[0,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVAnalytical, minStressDiv, maxStressDiv, xMin, xMax, yMin, yMax, \
                   False, False, r'$(\nabla \cdot \sigma)_{v^\prime}$ Analytical', r'(d)$\times20$', True)

    plot_subfigure(axes[1,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Wachs. ($u^\prime$ direction)', '(e)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVWachDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Wachs. ($v^\prime$ direction)', '(f)', False)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUPWLDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'PWL ($u^\prime$ direction)', '(g)', False)
    plot_subfigure(axes[1,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVPWLDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'PWL ($v^\prime$ direction)', '(h)', True)

    plot_subfigure(axes[2,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUWeakWachDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak Wachs. ($u^\prime$ direction)', '(i)', False)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVWeakWachDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak Wachs. ($v^\prime$ direction)', '(j)', False)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUWeakPWLDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak PWL ($u^\prime$ direction)', '(k)', False)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVWeakPWLDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak PWL ($v^\prime$ direction)', '(l)', True)

    plot_subfigure(axes[3,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak ($u^\prime$ direction)', '(m)', False)
    plot_subfigure(axes[3,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, latVertex, stressDivergenceVWeakDiff, minStressDivDiff, maxStressDivDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'Weak ($v^\prime$ direction)', '(n)', True)
    axes[3,2].axis('off')
    axes[3,3].axis('off')



    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_stress_divergence_map.png",dpi=400)
    #plt.savefig("strain_stress_divergence_map_3.png", bbox_inches="tight",dpi=2000)

    plt.clf()
    plt.cla()
    plt.close()



#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_stress_divergence_map()
