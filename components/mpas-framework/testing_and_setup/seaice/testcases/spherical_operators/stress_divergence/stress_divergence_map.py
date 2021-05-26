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
            for iCellOnVertex in range(0,vertexDegree):

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

def stress_divergence_map():

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

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    print("Stress divergence: ",
          np.amin(stressDivergenceUAnalytical), np.amax(stressDivergenceUAnalytical),
          np.amin(stressDivergenceVAnalytical), np.amax(stressDivergenceVAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_wachspress_40962/output.2000.nc","r")

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("Wachs: ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./output_pwl_40962/output.2000.nc","r")

    stressDivergenceUPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLDiff = (stressDivergenceUPWL - stressDivergenceUAnalytical)
    stressDivergenceVPWLDiff = (stressDivergenceVPWL - stressDivergenceVAnalytical)

    print("PWL:   ",
          np.amin(stressDivergenceUPWLDiff), np.amax(stressDivergenceUPWLDiff),
          np.amin(stressDivergenceVPWLDiff), np.amax(stressDivergenceVPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./output_weak_40962/output.2000.nc","r")

    stressDivergenceUWeak = fileWeak.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeak = fileWeak.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakDiff = (stressDivergenceUWeak - stressDivergenceUAnalytical)
    stressDivergenceVWeakDiff = (stressDivergenceVWeak - stressDivergenceVAnalytical)

    print("Weak:  ",
          np.amin(stressDivergenceUWeakDiff), np.amax(stressDivergenceUWeakDiff),
          np.amin(stressDivergenceVWeakDiff), np.amax(stressDivergenceVWeakDiff))

    fileWeak.close()


    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minVelocity = -1.0
    maxVelocity =  1.0

    minStressDiv = -20.0
    maxStressDiv =  20.0

    minStressDivDiff = -1.0
    maxStressDivDiff =  1.0


    fig, axes = plt.subplots(4, 4)

    fig.set_size_inches(7, 6.75)

    plot_subfigure(axes[0,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, minVelocity, maxVelocity, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Analytical', r'(b)$\times20$', False)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, minVelocity, maxVelocity, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(c)', False)
    plot_subfigure(axes[0,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{v^\prime}$ Analytical', r'(d)$\times20$', True)

    plot_subfigure(axes[1,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Wachs. ($u^\prime$ direction)', '(e)', False)
    plot_subfigure(axes[2,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'PWL ($u^\prime$ direction)', '(i)', False)
    plot_subfigure(axes[3,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Weak ($u^\prime$ direction)', '(m)', False)

    plot_subfigure(axes[1,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'Wachs. ($u^\prime$ direction)', '(f)', False)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUPWLDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'PWL ($u^\prime$ direction)', '(j)', False)
    plot_subfigure(axes[3,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'Weak ($u^\prime$ direction)', '(n)', False)

    plot_subfigure(axes[1,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWachDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Wachs. ($v^\prime$ direction)', '(g)', False)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVPWLDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'PWL ($v^\prime$ direction)', '(k)', False)
    plot_subfigure(axes[3,2], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWeakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'Weak ($v^\prime$ direction)', '(o)', False)

    plot_subfigure(axes[1,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWachDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'Wachs. ($v^\prime$ direction)', '(h)', True)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVPWLDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'PWL ($v^\prime$ direction)', '(l)', True)
    plot_subfigure(axes[3,3], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWeakDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'Weak ($v^\prime$ direction)', '(p)', True)

    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_divergence_map.png",dpi=400)
    #plt.savefig("stress_divergence_map_3.png", bbox_inches="tight",dpi=2000)

    plt.clf()
    plt.cla()
    plt.close()


    # histograms

    stressDivergenceUWachDiffHist = []
    stressDivergenceUPWLDiffHist  = []
    stressDivergenceUWeakDiffHist = []

    maxValue = 0.0
    for iVertex in range(0,nVertices):

        if (latVertex[iVertex] > 20.0*degreesToRadians):

            stressDivergenceUWachDiffHist.append(math.fabs(stressDivergenceUWachDiff[iVertex]))
            stressDivergenceUPWLDiffHist.append(math.fabs(stressDivergenceUPWLDiff[iVertex]))
            stressDivergenceUWeakDiffHist.append(math.fabs(stressDivergenceUWeakDiff[iVertex]))

            maxValue = max(math.fabs(stressDivergenceUWachDiff[iVertex]),maxValue)
            maxValue = max(math.fabs(stressDivergenceUPWLDiff[iVertex]),maxValue)


    mpl.rc('text', usetex=True)
    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rcParams['axes.linewidth'] = 0.5

    plt.figure(figsize=(3.74016, 3))

    #plt.hist(stressDivergenceUWachDiffHist, 50, normed=0, range=[0.0,0.08], histtype='step', lw=1, color='blue',  label='Wachspress')
    #plt.hist(stressDivergenceUPWLDiffHist,  50, normed=0, range=[0.0,0.08], histtype='step', lw=1, color='red',   label='PWL')
    #plt.hist(stressDivergenceUWeakDiffHist, 50, normed=0, range=[0.0,0.08], histtype='step', lw=1, color='green', label='Weak')
    plt.hist(stressDivergenceUWachDiffHist, 50, range=[0.0,0.08], histtype='step', lw=1, color='blue',  label='Wachspress')
    plt.hist(stressDivergenceUPWLDiffHist,  50, range=[0.0,0.08], histtype='step', lw=1, color='red',   label='PWL')
    plt.hist(stressDivergenceUWeakDiffHist, 50, range=[0.0,0.08], histtype='step', lw=1, color='green', label='Weak')

    plt.yscale('log', nonpositive='clip')

    plt.xlabel("Error")
    plt.ylabel("Frequency")

    plt.legend(["Wachspress","PWL","Weak"], frameon=False, fontsize=8)

    plt.xlim([0,0.08])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    plt.savefig("stress_divergence_hist.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_divergence_map()
