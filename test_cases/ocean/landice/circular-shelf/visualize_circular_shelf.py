#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from optparse import OptionParser
import matplotlib.pyplot as plt
# from matplotlib.contour import QuadContourSet


parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize", metavar="FILE")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")
# parser.add_option("-v", "--var", dest="variable", help="variable to visualize", metavar="VAR")
# parser.add_option("--max", dest="maximum", help="maximum for color bar", metavar="MAX")
# parser.add_option("--min", dest="minimum", help="minimum for color bar", metavar="MIN")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

time_slice = 0


secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!


f = NetCDFFile(options.filename,'r')

times = f.variables['xtime']
thickness = f.variables['thickness']
dcedge = f.variables['dcEdge']
try:
   bedTopography = f.variables['bedTopography']  # not needed
except:
   print "bedTopography not in file.  Continuing without it."
xCell = f.variables['xCell'][:]/1000.0
yCell = f.variables['yCell'][:]/1000.0
xEdge = f.variables['xEdge'][:]/1000.0
yEdge = f.variables['yEdge'][:]/1000.0
angleEdge = f.variables['angleEdge']
lowerSurface = f.variables['lowerSurface']
upperSurface = f.variables['upperSurface']
normalVelocity = f.variables['normalVelocity']
uReconstructX = f.variables['uReconstructX']
uReconstructY = f.variables['uReconstructY']
layerInterfaceSigma = f.variables['layerInterfaceSigma'][:]

vert_levs = len(f.dimensions['nVertLevels'])

time_length = times.shape[0]


velnorm = (uReconstructX[:]**2 + uReconstructY[:]**2)**0.5 * secInYr
print "Maximum velocity (m/yr) at cell centers in domain:", velnorm.max()


##################
# FIGURE: Map view surface and bed velocity
##################
fig = plt.figure(1)
ax = fig.add_subplot(121, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, velnorm[time_slice,:,0], marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xCell[:], yCell[:], uReconstructX[time_slice,:, 0] * secInYr, uReconstructY[time_slice,:, 0] * secInYr )
plt.title('surface speed (m/yr)' )
plt.draw()
ax = fig.add_subplot(122, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, velnorm[time_slice,:,-1], marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xCell[:], yCell[:], uReconstructX[time_slice,:, -1] * secInYr, uReconstructY[time_slice,:, -1] * secInYr )
plt.title('basal speed (m/yr)' )
plt.draw()
if options.saveimages:
        plt.savefig('circ_shelf_velos.png')



##################
# FIGURE: Cross-section of surface velocity through centerline
##################
fig = plt.figure(2)
# find cells on a longitudinal cross-section at y=0
# Note this may not work on a variable resolution mesh
# Note this assumes the setup script put the center of a cell at the center of the mesh at 0,0
indXsect = np.where(yCell == 0.0)[0]
indXsectIce = np.where(np.logical_and(yCell == 0.0, thickness[time_slice,:]>0.0))[0]

# contour speed across the cross-section
plt.contourf(np.tile(xCell[indXsectIce], (vert_levs+1,1)).transpose(), 
   (np.tile(thickness[time_slice, indXsectIce], (vert_levs+1,1))* np.tile(layerInterfaceSigma, (len(indXsectIce),1)).transpose() + lowerSurface[time_slice, indXsectIce]).transpose(), 
   velnorm[time_slice,indXsectIce,:], 100)
# plot x's at the velocity data locations
plt.plot(np.tile(xCell[indXsectIce], (vert_levs+1,1)).transpose(), 
   (np.tile(thickness[time_slice, indXsectIce], (vert_levs+1,1))* np.tile(layerInterfaceSigma, (len(indXsectIce),1)).transpose() + lowerSurface[time_slice, indXsectIce]).transpose(),
    'kx')#, label='velocity points')
cbar=plt.colorbar()
cbar.set_label('speed (m/yr)', rotation=270)

plt.plot(xCell[indXsectIce], upperSurface[time_slice, indXsectIce], 'ro-', label="Upper surface")
plt.plot(xCell[indXsectIce], lowerSurface[time_slice, indXsectIce], 'bo-', label="Lower surface")
try:
   plt.plot(xCell[indXsect], bedTopography[time_slice, indXsect], 'go-', label="Bed topography")
except:
   print "Skipping plotting of bedTopography."
plt.plot(xCell[indXsect], xCell[indXsect] * 0.0, ':k', label="sea level")
plt.legend(loc='best')
plt.title('cross-section at y=0' )
plt.draw()
if options.saveimages:
        plt.savefig('circ_shelf_xsect.png')


##################
# FIGURE: scatter plot of upper and lower surfaces
##################
#fig = plt.figure(3)
#ax = fig.add_subplot(121, aspect='equal')
#plt.scatter(xCell[:], yCell[:], 80, upperSurface[time_slice,:], marker='h', edgecolors='none')
#plt.colorbar()
#plt.title('upper surface (m)' )
#plt.draw()
#ax = fig.add_subplot(122, aspect='equal')
#plt.scatter(xCell[:], yCell[:], 80, lowerSurface[time_slice,:], marker='h', edgecolors='none')
#plt.colorbar()
#plt.title('lower surface (m)' )
#plt.draw()
#if options.saveimages:
#        plt.savefig('circ_shelf_surfaces.png')


if options.hidefigs:
     print "Plot display disabled with -n argument."
else:     
     plt.show()

f.close()

