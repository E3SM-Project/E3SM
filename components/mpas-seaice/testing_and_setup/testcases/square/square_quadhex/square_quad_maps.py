from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
#from matplotlib.mlab import griddata
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import griddata

#---------------------------------------------------------------

def covert_cice_array_to_mpas_array(nVertices,ciceArray):

    minArray = 1e30
    maxArray = -1e30

    mpasArray = np.zeros(nVertices)

    for i in range(0,81):
        for j in range(0,81):

            i2 = i + 1
            j2 = j + 1

            ij = i + j * 81

            mpasArray[ij] = ciceArray[j2,i2]

            if (mpasArray[ij] < minArray): minArray = mpasArray[ij]
            if (mpasArray[ij] > maxArray): maxArray = mpasArray[ij]

    return mpasArray, minArray, maxArray

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#---------------------------------------------------------------

def get_mpas_patch_collection(nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, mpasArray, cmap, vmin, vmax, limitPatches=False, rasterizePatches=False):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    nPatches = 0
    for iVertex in range(0,nVertices):

        if (interiorVertex[iVertex] == 1 and (math.fabs(mpasArray[iVertex]) > 1e-3 or not limitPatches)):

            polygonVertices = []

            for iVertexDegree in range(0,vertexDegree):

                iCell = cellsOnVertex[iVertex,iVertexDegree] - 1

                polygonVertices.append((xCell[iCell],yCell[iCell]))

            polygon = Polygon(polygonVertices)
            patches.append(polygon)

            colours.append(mpasArray[iVertex])

            minval = min(minval,mpasArray[iVertex])
            maxval = max(maxval,mpasArray[iVertex])

            nPatches = nPatches + 1

    print("nPatches: ", nPatches, limitPatches)

    patchCollection = PatchCollection(patches, cmap=cmap, rasterized=rasterizePatches)
    patchCollection.set_array(np.array(colours))
    patchCollection.set_linewidth(0)
    patchCollection.set_clim(vmin=vmin, vmax=vmax)

    return patchCollection, minval, maxval

#---------------------------------------------------------------

def set_up_axes(axes, title=None, subfigureLabel=None):

    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    #axes.set_ylim([minY+0.75*(maxY-minY), minY+1.0*(maxY-minY)])
    #axes.set_xlim([minX+0.75*(maxX-minX), minX+1.0*(maxX-minX)])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

    if (title != None):
        #axes.set_title(title, y=1.06, fontsize=8)
        axes.set_title(subfigureLabel+" "+title, loc="left")

    #if (subfigureLabel != None):
    #    axes.text(0.12, 0.9, subfigureLabel, verticalalignment='bottom', horizontalalignment='right',transform=axes.transAxes, fontsize=8)

#---------------------------------------------------------------

def plot_subfigure_patch(axes, nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, array, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, vmin=None, vmax=None, colourbarTickLocs=None, limitPatches=False, rasterizePatches=False, colorbar=True):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.jet
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    # patch plot
    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, array, colourMap, vmin, vmax, limitPatches=limitPatches, rasterizePatches=rasterizePatches)
    axes.add_collection(patchCollection)

    set_up_axes(axes, title, subfigureLabel)

    if (colorbar):
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        if (colourbarTickLocs != None):
            cb = plt.colorbar(patchCollection,cax=cax,ticks=colourbarTickLocs)
        else:
            cb = plt.colorbar(patchCollection,cax=cax)

        if (sciNote):
            cb.formatter.set_powerlimits((0, 0))
            cb.update_ticks()
        #if (diffPlot):
        #    cmLimit = max(math.fabs(minArray),math.fabs(maxArray))
        #   cb.set_clim(-cmLimit,cmLimit)
        #if (vmin != None):
        #    cb.set_clim(vmin,vmax)

    return patchCollection

#---------------------------------------------------------------

def get_contour_colours(nColours):

    contourColors = []

    for iColour in range(0,nColours):

        rColour = float(iColour) / float(nColours-1)
        iColour = 1.0 - rColour

        colour1 = [256, 0, 0]
        colour2 = [0, 0, 256]

        rColour1 = [float(colour1[0])/256.0,float(colour1[1])/256.0,float(colour1[2])/256.0]
        rColour2 = [float(colour2[0])/256.0,float(colour2[1])/256.0,float(colour2[2])/256.0]

        r = rColour * rColour1[0] + iColour * rColour2[0]
        g = rColour * rColour1[1] + iColour * rColour2[1]
        b = rColour * rColour1[2] + iColour * rColour2[2]

        contourColors.append((r,g,b))

    return tuple(contourColors)

#---------------------------------------------------------------

def plot_subfigure_contour(axes, nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, array, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, vmin=None, vmax=None, contourLevels=None, clabelLocations=None):

    # contour plot
    x = []
    y = []
    z = []
    points = []
    for iVertex in range(0,nVertices):
        if (interiorVertex[iVertex] == 1):

            xAvg = 0.0
            yAvg = 0.0

            for iVertexDegree in range(0,vertexDegree):

                iCell = cellsOnVertex[iVertex,iVertexDegree] - 1

                xAvg = xAvg + xCell[iCell]
                yAvg = yAvg + yCell[iCell]

            xAvg = xAvg / vertexDegree
            yAvg = yAvg / vertexDegree

            x.append(xAvg)
            y.append(yAvg)
            points.append([xAvg,yAvg])

            z.append(array[iVertex])

    xi = np.linspace(minX,maxX, 100)
    yi = np.linspace(minY,maxY, 100)
    #zi = griddata(x, y, z, xi, yi, interp='linear')

    pointsOut = []
    for i in range(0,100):
        for j in range(0,100):
            pointsOut.append([xi[i],yi[j]])
    zi = griddata(points, z, pointsOut, 'linear')
    zii = np.zeros((100,100))
    for i in range(0,100):
        for j in range(0,100):
            ij = j*100 + i
            zii[i,j] = zi[ij]


    contourPlot = axes.contour(xi, yi, zii, 10, linewidths=0.3, levels=contourLevels, colors='k', linestyles="solid")#get_contour_colours(len(contourLevels)))
    axes.clabel(contourPlot, inline=1, fontsize=5, fmt=FormatStrFormatter('%g'), manual=clabelLocations, inline_spacing=1)

    set_up_axes(axes, title, subfigureLabel)

    #divider = make_axes_locatable(axes)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #cax.set_axis_off()

    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    #axes.set_ylim([minY+0.75*(maxY-minY), minY+1.0*(maxY-minY)])
    #axes.set_xlim([minX+0.75*(maxX-minX), minX+1.0*(maxX-minX)])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

#---------------------------------------------------------------

def ml(x, y):

    xx = minX + x * (maxX - minX)
    yy = minY + y * (maxY - minY)

    return (xx,yy)

#---------------------------------------------------------------

print(mpl.__version__)

mpasResultsDirname = "./"
#ciceResultsDirname = "/home/akt/Work/MPAS-Seaice/papers/MPAS-seaice_intro_2/GMD_v1/figures/square/results/cice_internal_stress"
ciceResultsDirname = "/home/akt/Work/CICE/square/cice_2/CICE-svn-trunk-cice-5.1.2/rundirs/rundir/history/"
#ciceResultsDirname = "/home/akt/Work/CICE/square/CICE/rundir/history/"

#fileCICEstd = "iceh_inst.1997-01-01-03600.nc"
#fileCICEstd = "iceh_inst.1997-01-01-07200.nc"
#fileCICEstd = "iceh_inst.1997-01-01-10800.nc"
fileCICEstd = "iceh_inst.1997-01-01-14400.nc"

fileMPASstd = "output.2000.nc"
iTimeMPAS = 3

filenameCICE      = ciceResultsDirname + "/" + fileCICEstd
filenameMPASWachs = mpasResultsDirname + "output_quad_wachspress_0080x0080_120" + "/" + fileMPASstd
filenameMPASPWL   = mpasResultsDirname + "output_quad_pwl_0080x0080_120"        + "/" + fileMPASstd
filenameMPASWeak  = mpasResultsDirname + "output_quad_weak_0080x0080_120"       + "/" + fileMPASstd

fileCICE      = Dataset(filenameCICE,      "r")
fileMPASWachs = Dataset(filenameMPASWachs, "r")
fileMPASPWL   = Dataset(filenameMPASPWL,   "r")
fileMPASWeak  = Dataset(filenameMPASWeak,  "r")

# load CICE
uCICEin   = fileCICE.variables["uvel"][0,:]
divCICEin = fileCICE.variables["strintx"][0,:]

# Wachspress
uMPASWachs   = fileMPASWachs.variables["uVelocity"][iTimeMPAS,:]
divMPASWachs = fileMPASWachs.variables["stressDivergenceU"][iTimeMPAS,:]

# load MPAS geometry
nVertices    = len(fileMPASWachs.dimensions["nVertices"])
vertexDegree = len(fileMPASWachs.dimensions["vertexDegree"])

interiorVertex = fileMPASWachs.variables["interiorVertex"][iTimeMPAS,:]
cellsOnVertex  = fileMPASWachs.variables["cellsOnVertex"][:,:]
xCell          = fileMPASWachs.variables["xCell"][:]
yCell          = fileMPASWachs.variables["yCell"][:]
maxX = np.amax(xCell)
minX = np.amin(xCell)
maxY = np.amax(yCell)
minY = np.amin(yCell)

# PWL
uMPASPWL   = fileMPASPWL.variables["uVelocity"][iTimeMPAS,:]
divMPASPWL = fileMPASPWL.variables["stressDivergenceU"][iTimeMPAS,:]

# Weak
uMPASWeak   = fileMPASWeak.variables["uVelocity"][iTimeMPAS,:]
divMPASWeak = fileMPASWeak.variables["stressDivergenceU"][iTimeMPAS,:]

# convert cice arrays
uCICE,   uCICEMin,   uCICEMax   = covert_cice_array_to_mpas_array(nVertices,uCICEin)
divCICE, divCICEMin, divCICEMax = covert_cice_array_to_mpas_array(nVertices,divCICEin)

print("CICE u:       ", np.nanmin(uCICE),    np.nanmax(uCICE))
print("MPAS Wachs u: ", np.amin(uMPASWachs), np.amax(uMPASWachs))
print("MPAS PWL u:   ", np.amin(uMPASPWL),   np.amax(uMPASPWL))
print("MPAS Weak u:  ", np.amin(uMPASWeak),  np.amax(uMPASWeak))


# difference arrays
uDiffWachs = uMPASWachs - uCICE
uDiffPWL   = uMPASPWL   - uCICE
uDiffWeak  = uMPASWeak  - uCICE

divDiffWachs = divMPASWachs - divCICE
divDiffPWL   = divMPASPWL   - divCICE
divDiffWeak  = divMPASWeak  - divCICE

fileCICE.close()
fileMPASWachs.close()
fileMPASPWL.close()
fileMPASWeak.close()


print("nVertices:    ", nVertices)
print("vertexDegree: ", vertexDegree)




uCICEMin = np.nanmin(uCICE)
uCICEMax = np.nanmax(uCICE)

divCICEMin = np.nanmin(divCICE)
divCICEMax = np.nanmax(divCICE)

# rescale
uDiffWachs = uDiffWachs / uCICEMax
uDiffPWL   = uDiffPWL   / uCICEMax
uDiffWeak  = uDiffWeak  / uCICEMax



uDiffWachsMin = np.nanmin(uDiffWachs)
uDiffWachsMax = np.nanmax(uDiffWachs)

uDiffPWLMin = np.nanmin(uDiffPWL)
uDiffPWLMax = np.nanmax(uDiffPWL)

uDiffWeakMin = np.nanmin(uDiffWeak)
uDiffWeakMax = np.nanmax(uDiffWeak)

divDiffWachsMin = np.nanmin(divDiffWachs)
divDiffWachsMax = np.nanmax(divDiffWachs)

divDiffPWLMin = np.nanmin(divDiffPWL)
divDiffPWLMax = np.nanmax(divDiffPWL)

divDiffWeakMin = np.nanmin(divDiffWeak)
divDiffWeakMax = np.nanmax(divDiffWeak)




uMin = uCICEMin
uMax = uCICEMax

divMin = divCICEMin
divMax = divCICEMax

uDiffMin = min(uDiffWachsMin, uDiffPWLMin, uDiffWeakMin)
uDiffMax = max(uDiffWachsMax, uDiffPWLMax, uDiffWeakMax)

divDiffMin = min(divDiffWachsMin, divDiffPWLMin, divDiffWeakMin)
divDiffMax = max(divDiffWachsMax, divDiffPWLMax, divDiffWeakMax)

print(uMin, uMax)
print(divMin, divMax)

print(uDiffMin, uDiffMax)
print(divDiffMin, divDiffMax)

# fix for final plot
uMin = min(-0.0151097934221, uMin)
uMax = max(0.182248369745, uMax)

divMin = min(-0.196414155474, divMin)
divMax = max(0.0444513971139, divMax)

#uMin = min(-0.02, uMin)
#uMax = max(0.22, uMax)

#divMin = min(-0.196414155474, divMin)
#divMax = max(0.0444513971139, divMax)

uDiffMin = -0.02
uDiffMax = 0.02

divDiffMin = -0.08
divDiffMax = 0.08

uMin = -0.00593214085963
uMax = 0.205460178207

divMin = -0.2#-0.191937435754
divMax = 0.05#0.0474654996112


# plot
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


lx = 15
ly = 7.8
fig, axes = plt.subplots(2,4,figsize=(lx*cm,ly*cm))
rasterizePatches = True



# AGU sizing
plot_subfigure_contour(axes[0,0], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, uCICE       , \
                       False, \
                       False, \
                       r'$u$ CICE', \
                       "(a)", \
                       uMin, uMax, \
                       contourLevels=[0.0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2], \
                       clabelLocations = [ml(0.95,0.2),ml(0.45,0.1),ml(0.1,0.2),ml(0.55,0.3),ml(0.45,0.5),ml(0.55,0.7),ml(0.5,0.9)])

plot_subfigure_contour(axes[1,0], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, divCICE     , \
                       False, \
                       False, \
                       r'$(\nabla \cdot \sigma)_{u}$ CICE', \
                       "(e)", \
                       divMin, divMax, \
                       contourLevels=[-0.1, -0.075, -0.05, -0.02, -0.01, -0.005, 0])

im0 = plot_subfigure_patch(axes[0,1], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, uDiffWachs  , \
                           True, \
                           True, \
                           r'$\Delta u/u_{max}$ Wachspress', \
                           "(b)", \
                           uDiffMin, uDiffMax, \
                           limitPatches=True, \
                           rasterizePatches=rasterizePatches, \
                           colorbar=False)

im1 = plot_subfigure_patch(axes[1,1], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, divDiffWachs, \
                           True, \
                           True, \
                           r'$\Delta (\nabla \cdot \sigma)_{u}$ Wachspress', \
                           "(f)", \
                           divDiffMin, divDiffMax,\
                           limitPatches=True, \
                           rasterizePatches=rasterizePatches, \
                           colorbar=False)

plot_subfigure_patch(axes[0,2], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, uDiffPWL    , \
                     True, \
                     True, \
                     r'$\Delta u/u_{max}$ PWL', \
                     "(c)", \
                     uDiffMin, uDiffMax, \
                     limitPatches=True, \
                     rasterizePatches=rasterizePatches, \
                     colorbar=False)

plot_subfigure_patch(axes[1,2], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, divDiffPWL  , \
                     True, \
                     True, \
                     r'$\Delta (\nabla \cdot \sigma)_{u}$ PWL', \
                     "(g)", \
                     divDiffMin, divDiffMax, \
                     limitPatches=True, \
                     rasterizePatches=rasterizePatches, \
                     colorbar=False)

plot_subfigure_patch(axes[0,3], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, uDiffWeak   , \
                     True, \
                     True, \
                     r'$\Delta u/u_{max}$ Weak', \
                     "(d)", \
                     uDiffMin, uDiffMax, \
                     limitPatches=True, \
                     rasterizePatches=rasterizePatches, \
                     colorbar=False)

plot_subfigure_patch(axes[1,3], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, divDiffWeak , \
                     True, \
                     True, \
                     r'$\Delta (\nabla \cdot \sigma)_{u}$ Weak', \
                     "(h)", \
                     divDiffMin, divDiffMax, \
                     limitPatches=True, \
                     rasterizePatches=rasterizePatches, \
                     colorbar=False)

mainWidthThreePanel = axes[0,2].get_position().x1 * lx
print("Three panel width: ", mainWidthThreePanel)


mainWidth = 13.5
leftMargin = 0.3

left = leftMargin / lx
right = mainWidth / lx
dw = 0.75
dh = 0.75
plt.subplots_adjust(left=left,
                    bottom=0.01,
                    right=right,
                    top=0.95,
                    wspace=dw/lx,
                    hspace=dh/ly)




cbwidth = 0.25
cbdx = cbwidth / lx
cbx0 = right + cbdx


pos0 = axes[0,0].get_position()
cb_ax = fig.add_axes([cbx0, pos0.y0, cbdx, pos0.y1-pos0.y0])
cb = fig.colorbar(im0, cax=cb_ax)

pos1 = axes[1,0].get_position()
cb_ax = fig.add_axes([cbx0, pos1.y0, cbdx, pos1.y1-pos1.y0])
cb = fig.colorbar(im1, cax=cb_ax)


plt.savefig("square_quad_maps.eps",dpi=300)
plt.savefig("square_quad_maps.png",dpi=300)


# analysis
fig, axes = plt.subplots(2,3,figsize=(11.25*cm,7.4*cm))

print(np.nanmin(uCICE),np.nanmax(uCICE))
print(np.nanmin(uMPASWachs),np.nanmax(uMPASWachs))


# AGU sizing
plot_subfigure_patch(axes[0,0], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, uCICE       , \
               False, \
               False, \
               r'$u$ CICE', \
               "(a)", \
               min(np.nanmin(uCICE),np.nanmin(uMPASWachs)), max(np.nanmax(uCICE),np.nanmax(uMPASWachs)), \
               limitPatches=False, \
                     colorbar=True)

plot_subfigure_patch(axes[1,0], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, divCICE     , \
               False, \
               False, \
               r'$(\nabla \cdot \sigma)_{u}$ CICE', \
               "(d)", \
               min(np.nanmin(divCICE),np.nanmin(divMPASWachs)), max(np.nanmax(divCICE),np.nanmax(divMPASWachs)), \
               limitPatches=False, \
                     colorbar=True)

plot_subfigure_patch(axes[0,1], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, uMPASWachs  , \
               False, \
               False, \
               r'$u$ Wachspress', \
               "(b)", \
               min(np.nanmin(uCICE),np.nanmin(uMPASWachs)), max(np.nanmax(uCICE),np.nanmax(uMPASWachs)), \
               limitPatches=False, \
                     colorbar=True)

plot_subfigure_patch(axes[1,1], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, divMPASWachs, \
               False, \
               False, \
               r'$\Delta (\nabla \cdot \sigma)_{u}$ Wachspress', \
               "(e)", \
               min(np.nanmin(divCICE),np.nanmin(divMPASWachs)), max(np.nanmax(divCICE),np.nanmax(divMPASWachs)), \
               limitPatches=False, \
                     colorbar=True)

plot_subfigure_patch(axes[0,2], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, uDiffWachs    , \
               True, \
               True, \
               r'$\Delta u$ Wachspress', \
               "(c)", \
               np.nanmin(uDiffWachs), np.nanmax(uDiffWachs), \
               limitPatches=False, \
                     colorbar=True)

plot_subfigure_patch(axes[1,2], nVertices, vertexDegree, interiorVertex, cellsOnVertex, xCell, yCell, divDiffWachs  , \
               True, \
               True, \
               r'$\Delta (\nabla \cdot \sigma)_{u}$ Wachspress', \
               "(f)", \
               np.nanmin(divDiffWachs), np.nanmax(divDiffWachs), \
               limitPatches=False, \
                     colorbar=True)


plt.savefig("square_quad_comparison.eps",dpi=300)
plt.savefig("square_quad_comparison.png",dpi=300)
