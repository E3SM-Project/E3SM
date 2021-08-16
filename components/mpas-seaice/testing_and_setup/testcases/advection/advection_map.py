from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math, sys
from scipy.interpolate import griddata
from matplotlib.ticker import FormatStrFormatter

#---------------------------------------------------------------

def get_mpas_patch_collection(nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, mpasArray, cmap):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iCell in range(0,nCells):

        polygonVertices = []

        lUse = True

        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

            iVertex = verticesOnCell[iCell,iVertexOnCell] - 1

            polygonVertices.append((xVertex[iVertex],zVertex[iVertex]))

            if (yVertex[iVertex] < 0.0):
                lUse = False

        if (lUse):

            polygon = Polygon(polygonVertices)
            patches.append(polygon)

            colours.append(mpasArray[iCell])

            minval = min(minval,mpasArray[iCell])
            maxval = max(maxval,mpasArray[iCell])

    patchCollection = PatchCollection(patches, cmap=cmap, rasterized=False)
    patchCollection.set_array(np.array(colours))
    patchCollection.set_linewidth(0)
    patchCollection.set_clim(vmin=0.0,vmax=1.0)

    return patchCollection, minval, maxval

#---------------------------------------------------------------

def plot_subfigure_patch(axes,
                         nCells,
                         nEdgesOnCell,
                         verticesOnCell,
                         xVertex,
                         yVertex,
                         zVertex,
                         array,
                         minX,
                         maxX,
                         minY,
                         maxY,
                         sciNote=False,
                         diffPlot=False,
                         title=None,
                         cbTitle=None,
                         subfigureLabel=None):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.plasma
    else:
        colourMap = mpl.cm.RdBu

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, array, colourMap)
    axes.add_collection(patchCollection)
    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

    if (title != None):
        axes.set_title(title)

    if (subfigureLabel != None):
        axes.text(0.12, 0.9, subfigureLabel, verticalalignment='bottom', horizontalalignment='right',transform=axes.transAxes, fontsize=8)

    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(patchCollection,cax=cax)
    #if (cbTitle != None):
    #    cb.ax.set_ylabel(cbTitle)
    if (sciNote):
        cb.formatter.set_powerlimits((0, 0))
        cb.update_ticks()
    if (diffPlot):
        cmLimit = max(math.fabs(minArray),math.fabs(maxArray))
        cb.set_clim(-cmLimit,cmLimit)

#---------------------------------------------------------------

def plot_subfigure_contour(axes,
                           nCells,
                           xCell,
                           yCell,
                           zCell,
                           array,
                           minX,
                           maxX,
                           minY,
                           maxY,
                           sciNote=False,
                           diffPlot=False,
                           title=None,
                           cbTitle=None,
                           subfigureLabel=None,
                           contourLevels=None,
                           addCLabels=False):

    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

    if (title != None):
        axes.set_title(title)

    if (subfigureLabel != None):
        #axes.text(0.12, 0.9, subfigureLabel, verticalalignment='bottom', horizontalalignment='right',transform=axes.transAxes, fontsize=8)
        axes.set_title(subfigureLabel, loc='left')

    # contour plot
    x = []
    y = []
    z = []
    for iCell in range(0,nCells):

        if (yCell[iCell] > 0.0):

            x.append(xCell[iCell])
            y.append(zCell[iCell])
            z.append(array[iCell])

    xi = np.linspace(minX, maxX, 100)
    yi = np.linspace(minY, maxY, 100)
    X, Y = np.meshgrid(xi,yi)

    print("  Griddata...")
    zi = griddata((x,y), z, (X,Y), method='linear', fill_value=0)

    print("  Contour plot...")
    contourPlot = axes.contour(X, Y, zi, 10, linewidths=0.5, colors='k', levels=contourLevels)
    if (addCLabels):
        axes.clabel(contourPlot, inline=1, fontsize=5, inline_spacing=1, fmt=FormatStrFormatter('%g'))#, manual=clabelLocations, inline_spacing=1)

#---------------------------------------------------------------

def advection_map():

    iTime = -1

    #res = "2562"
    #res = "10242"
    res = "40962"
    #res = "163842"

    experiment1 = "cosine_bell"
    experiment2 = "slotted_cylinder"

    print("Load data...")

    # slotted cylinder IR
    filenameIR     = "./output_IR_%s_%s/output.2000.nc" %(experiment2,res)

    fileIR = Dataset(filenameIR, "r")

    nCells = len(fileIR.dimensions["nCells"])

    #sphereRadius = fileIR.attributes["sphere_radius"]
    sphereRadius = 6371229.0

    nEdgesOnCell   = fileIR.variables["nEdgesOnCell"][:]
    verticesOnCell = fileIR.variables["verticesOnCell"][:]
    xVertex        = fileIR.variables["xVertex"][:]
    yVertex        = fileIR.variables["yVertex"][:]
    zVertex        = fileIR.variables["zVertex"][:]
    xCell          = fileIR.variables["xCell"][:]
    yCell          = fileIR.variables["yCell"][:]
    zCell          = fileIR.variables["zCell"][:]

    iceAreaInitialSC   = fileIR.variables["iceAreaCategory"][0,:,0,0]
    iceVolumeInitialSC = fileIR.variables["iceVolumeCategory"][0,:,0,0]

    iceAreaIRSC   = fileIR.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolumeIRSC = fileIR.variables["iceVolumeCategory"][iTime,:,0,0]

    fileIR.close()

    # slotted cylinder upwind
    filenameUpwind = "./output_upwind_%s_%s/output.2000.nc" %(experiment2,res)

    fileUpwind = Dataset(filenameUpwind, "r")

    iceAreaUpwindSC   = fileIR.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolumeUpwindSC = fileIR.variables["iceVolumeCategory"][iTime,:,0,0]

    fileUpwind.close()

    # cosine bell IR
    filenameIR = "./output_IR_%s_%s/output.2000.nc" %(experiment1,res)

    fileIR = Dataset(filenameIR, "r")

    iceAreaInitialCB   = fileIR.variables["iceAreaCategory"][0,:,0,0]
    iceVolumeInitialCB = fileIR.variables["iceVolumeCategory"][0,:,0,0]

    iceAreaIRCB   = fileIR.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolumeIRCB = fileIR.variables["iceVolumeCategory"][iTime,:,0,0]

    fileIR.close()

    # cosine bell Upwind
    filenameUpwind = "./output_upwind_%s_%s/output.2000.nc" %(experiment1,res)

    fileUpwind = Dataset(filenameUpwind, "r")

    iceAreaInitialCB   = fileUpwind.variables["iceAreaCategory"][0,:,0,0]
    iceVolumeInitialCB = fileUpwind.variables["iceVolumeCategory"][0,:,0,0]

    iceAreaUpwindCB   = fileUpwind.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolumeUpwindCB = fileUpwind.variables["iceVolumeCategory"][iTime,:,0,0]

    fileUpwind.close()


    scaleFactor = 0.9
    #scaleFactor = 1.0

    minX = -sphereRadius * scaleFactor
    maxX =  sphereRadius * scaleFactor

    minY = -sphereRadius * scaleFactor
    maxY =  sphereRadius * scaleFactor

    # plotting

    cm = 1/2.54  # centimeters in inches
    plt.rcParams["font.family"] = "Times New Roman"
    #mpl.rc('text', usetex=True)
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

    #mpl.rcParams['axes.linewidth'] = 0.5


    print("Create plots...")

    fig, axes = plt.subplots(2, 3, figsize=(15*cm,11*cm))

    contour = True

    if (contour):

        plot_subfigure_contour(axes[0,0], nCells, xCell, yCell, zCell, iceAreaInitialCB, minX, maxX, minY, maxY, False, False, r'', 'm',
                               '(a) CB initial condition', contourLevels=[0.05,0.1,0.3,0.5,0.7,0.9,0.95])
        plot_subfigure_contour(axes[1,0], nCells, xCell, yCell, zCell, iceAreaInitialSC, minX, maxX, minY, maxY, False, False, r'', 'm',
                               '(d) SC initial condition', contourLevels=[0.05,0.1,0.3,0.5,0.7,0.9,0.95])
        plot_subfigure_contour(axes[0,1], nCells, xCell, yCell, zCell, iceAreaUpwindCB,  minX, maxX, minY, maxY, False, False, r'', 'm',
                               '(b) CB upwind', contourLevels=[0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2], addCLabels=True)
        plot_subfigure_contour(axes[1,1], nCells, xCell, yCell, zCell, iceAreaUpwindSC,  minX, maxX, minY, maxY, False, False, r'', 'm',
                               '(e) SC upwind', contourLevels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], addCLabels=True)
        plot_subfigure_contour(axes[0,2], nCells, xCell, yCell, zCell, iceAreaIRCB,      minX, maxX, minY, maxY, False, False, r'', 'm',
                               '(c) CB incremental remapping', contourLevels=[0.05,0.1,0.3,0.5,0.7,0.9,0.95])
        plot_subfigure_contour(axes[1,2], nCells, xCell, yCell, zCell, iceAreaIRSC,      minX, maxX, minY, maxY, False, False, r'', 'm',
                               '(f) SC incremental remapping', contourLevels=[0.05,0.1,0.3,0.5,0.7,0.9,0.95])

        #plt.tight_layout()
        plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
        plt.savefig("advection_map.png",dpi=300)
        plt.savefig("advection_map.eps")

    else:

        plot_subfigure_patch(axes[0,0], nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, iceAreaInitialCB, minX, maxX, minY, maxY, False, False, r'', 'm', '(a)')
        plot_subfigure_patch(axes[1,0], nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, iceAreaInitialSC, minX, maxX, minY, maxY, False, False, r'', 'm', '(d)')
        plot_subfigure_patch(axes[0,1], nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, iceAreaUpwindCB,  minX, maxX, minY, maxY, False, False, r'', 'm', '(b)')
        plot_subfigure_patch(axes[1,1], nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, iceAreaUpwindSC,  minX, maxX, minY, maxY, False, False, r'', 'm', '(e)')
        plot_subfigure_patch(axes[0,2], nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, iceAreaIRCB,      minX, maxX, minY, maxY, False, False, r'', 'm', '(c)')
        plot_subfigure_patch(axes[1,2], nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, iceAreaIRSC,      minX, maxX, minY, maxY, False, False, r'', 'm', '(f)')

        plt.savefig("advection_map.png",bbox_inches="tight",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    advection_map()
