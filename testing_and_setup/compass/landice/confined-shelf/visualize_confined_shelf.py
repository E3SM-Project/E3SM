#!/usr/bin/env python

# Visualize results of the confined shelf experiment
# See http://homepages.vub.ac.be/~phuybrec/eismint/shelf-descr.pdf for description.

import numpy as np
from netCDF4 import Dataset as NetCDFFile
from optparse import OptionParser
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
# from matplotlib.contour import QuadContourSet


parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize", metavar="FILE")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

time_slice = 0


secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!


def contourMPAS(field, contour_levs):
  """Contours irregular MPAS data on cells"""
  #-- Now let's grid your data.
  # First we'll make a regular grid to interpolate onto. 
  numcols, numrows = nCells**0.5, nCells**0.5  # may want to adjust the density of the regular grid
  xi = np.linspace(xCell.min(), xCell.max(), numcols)
  yi = np.linspace(xCell.min(), yCell.max(), numrows)
  xi, yi = np.meshgrid(xi, yi)
  #-- Interpolate at the points in xi, yi
  zi = griddata( (xCell, yCell), field, (xi, yi) )
  #-- Display the results
  im = plt.contour(xi, yi, zi, contour_levs)
  #plt.scatter(xCell, yCell, c=temperature[timelev,:,-1], s=100, vmin=zi.min(), vmax=zi.max())  # to see the raw data on top
  plt.colorbar(im)




f = NetCDFFile(options.filename,'r')

times = f.variables['xtime']
thickness = f.variables['thickness']
dcedge = f.variables['dcEdge']
bedTopography = f.variables['bedTopography']  # not needed
xCell = f.variables['xCell'][:]/1000.0
yCell = f.variables['yCell'][:]/1000.0
xEdge = f.variables['xEdge'][:]/1000.0
yEdge = f.variables['yEdge'][:]/1000.0
angleEdge = f.variables['angleEdge']
temperature = f.variables['temperature']
lowerSurface = f.variables['lowerSurface']
upperSurface = f.variables['upperSurface']
normalVelocity = f.variables['normalVelocity']
uReconstructX = f.variables['uReconstructX']
uReconstructY = f.variables['uReconstructY']
layerCenterSigma = f.variables['layerCenterSigma'][:]
nCells = len(f.dimensions['nCells'])

vert_levs = len(f.dimensions['nVertLevels'])

time_length = times.shape[0]


velnorm = (uReconstructX[:]**2 + uReconstructY[:]**2)**0.5
print "Maximum speed (m/yr) at cell centers in domain:", velnorm.max() * secInYr

var_slice = thickness[time_slice,:]
# var_slice = var_slice.reshape(time_length, ny, nx)

##################
# FIGURE: Map view surface and bed velocity as scatter plot
##################
fig = plt.figure(1)
ax = fig.add_subplot(121, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, velnorm[time_slice,:,0] * secInYr, marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xCell[:], yCell[:], uReconstructX[time_slice,:, 0] * secInYr, uReconstructY[time_slice,:, 0] * secInYr )
plt.title('surface speed (m/yr)' )
plt.draw()
ax = fig.add_subplot(122, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, velnorm[time_slice,:,-1] * secInYr, marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xCell[:], yCell[:], uReconstructX[time_slice,:, -1] * secInYr, uReconstructY[time_slice,:, -1] * secInYr )
plt.title('basal speed (m/yr)' )
plt.draw()
if options.saveimages:
        plt.savefig('conf_shelf_velos.png')



##################
# FIGURE: Map view surface and bed velocity as contour plot
##################
fig = plt.figure(2, facecolor='w')

fig.add_subplot(121)
contourMPAS(uReconstructX[time_slice,:,0]*secInYr, np.linspace(-200.0, 200.0, 400.0/20.0+1))
plt.axis('equal')
plt.title('Horizontal x-velocity (m/a)')
plt.xlim( (0.0, 200.0) ); plt.ylim( (0.0, 200.0) )
plt.xlabel('X position (km)')
plt.ylabel('Y position (km)')

fig.add_subplot(122)
contourMPAS(uReconstructY[time_slice,:,0]*secInYr, np.linspace(-1000.0, 0.0, 1000.0/100.0+1))
plt.axis('equal')
plt.title('Horizontal y-velocity (m/a)')
plt.xlim( (0.0, 200.0) ); plt.ylim( (0.0, 200.0) )
plt.xlabel('X position (km)')
plt.ylabel('Y position (km)')

plt.draw()
print "Compare Figure 2 to test 3 results at http://homepages.vub.ac.be/~phuybrec/eismint/shelf-descr.pdf"
if options.saveimages:
        plt.savefig('conf_shelf_velo_contours.png')

##################
# FIGURE: Cross-section of surface velocity through centerline
##################
#fig = plt.figure(3)
## find cells on a longitudinal cross-section at y=0
## Note this may not work on a variable resolution mesh
## Note this assumes the setup script put the center of a cell at the center of the mesh at 0,0
#indXsect = np.where(yCell == 0.0)[0]
#indXsectIce = np.where(np.logical_and(yCell == 0.0, thickness[time_slice,:]>0.0))[0]

#plt.contourf(np.tile(xCell[indXsectIce], (vert_levs,1)).transpose(), 
#   (np.tile(thickness[time_slice, indXsectIce], (vert_levs,1))* np.tile(layerCenterSigma, (len(indXsectIce),1)).transpose() + lowerSurface[time_slice, indXsectIce]).transpose(), 
#   velnorm[time_slice,indXsectIce,:], 100)
#plt.plot(np.tile(xCell[indXsectIce], (vert_levs,1)).transpose(), 
#   (np.tile(thickness[time_slice, indXsectIce], (vert_levs,1))* np.tile(layerCenterSigma, (len(indXsectIce),1)).transpose() + lowerSurface[time_slice, indXsectIce]).transpose(),
#    'kx')#, label='velocity points')
#cbar=plt.colorbar()
#cbar.set_label('speed (m/yr)', rotation=270)

#plt.plot(xCell[indXsectIce], upperSurface[time_slice, indXsectIce], 'ro-', label="Upper surface")
#plt.plot(xCell[indXsectIce], lowerSurface[time_slice, indXsectIce], 'bo-', label="Lower surface")
#plt.plot(xCell[indXsect], bedTopography[indXsect], 'go-', label="Bed topography")
#plt.plot(xCell[indXsect], xCell[indXsect] * 0.0, ':k', label="sea level")
#plt.legend(loc='best')
#plt.title('cross-section at y=0' )
#plt.draw()
#if options.saveimages:
#        plt.savefig('conf_shelf_xsect.png')


##################
# FIGURE: scatter plot of upper and lower surfaces
##################
#fig = plt.figure(4)
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
#        plt.savefig('conf_shelf_surfaces.png')



if options.hidefigs:
     print "Plot display disabled with -n argument."
else:     
     plt.show()

f.close()

