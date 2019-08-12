#!/usr/bin/env python
import numpy
import netCDF4
import os
import scipy.ndimage.filters as filters
from mpl_toolkits.basemap import interp


def smoothGeometry(field, filterSigma):

    smoothedField = filters.gaussian_filter(field, filterSigma,
                                            mode='constant', cval=0.)

    return smoothedField

def readVar(varName, defaultValue=0.0):
    field = defaultValue*numpy.ones((ny,nx),float)
    field[buffer:-buffer,buffer:-buffer] = numpy.array(inFile.variables[varName])[:,minIndex:]
    return field

def writeVar(outVarName, field, attrs):
    outVar = outFile.createVariable(outVarName, 'f8', ('Time', 'nCells'))
    outVar[:, :] = field
    outVar.setncatts(attrs)


inIceGeomFileName = '../init_step1/input_geometry.nc'
inMpasMeshFileName = 'init.nc'
outFileName = 'land_ice_forcing.nc'
filterSigma = 2.0
initialScale = 0.1

buffer = 1

xMax = 800e3 #km
x0 = 320e3 #km

iceDensity = 918.
oceanDensity = 1028.
g = 9.80665
  
inFile = netCDF4.Dataset(inIceGeomFileName, 'r')
x = numpy.array(inFile.variables['x'])[:]
y = numpy.array(inFile.variables['y'])[:]

deltaX = x[1]-x[0]
deltaY = y[1]-y[0]

minIndex = numpy.nonzero(x >= x0)[0][0]
#print minIndex, x[minIndex]

nx = len(x)-minIndex+2*buffer
ny = len(y)+2*buffer

outX = numpy.array(nx)
outY = numpy.array(ny)

outX = x[minIndex] + deltaX*(-buffer + numpy.arange(nx))
#print outX
outY = y[0] + deltaY*(-buffer + numpy.arange(ny))

surf = readVar('upperSurface')
draft = readVar('lowerSurface')
groundedMask = readVar('groundedMask', defaultValue=1.0)
inFile.close()

iceThickness = surf-draft

icePressure = iceDensity*g*iceThickness

smoothedIcePressure = smoothGeometry(icePressure, filterSigma)

oceanFraction = 1. - groundedMask
smoothedMask = smoothGeometry(oceanFraction, filterSigma)
smoothedDraft = smoothGeometry(oceanFraction*draft, filterSigma)
threshold = 0.01
mask = smoothedMask > threshold
smoothedDraft[mask] /= smoothedMask[mask]

# interpolate
mpasFile = netCDF4.Dataset(inMpasMeshFileName, 'r')
xCell = mpasFile.variables['xCell'][:]
yCell = mpasFile.variables['yCell'][:]
mpasFile.close()

nCells = len(xCell)
StrLen = 64
xtime = numpy.zeros((2), 'S64')

xtime[0] = "0001-01-01_00:00:00                                             "
xtime[1] = "0002-01-01_00:00:00                                             "

landIcePressure = numpy.zeros((2, nCells), float)
landIcePressure[1, :] = interp(smoothedIcePressure, outX, outY, xCell, yCell)
landIceFraction = numpy.zeros(landIcePressure.shape)
landIceDraft = numpy.zeros(landIcePressure.shape)
landIceDraft[1, :] = interp(smoothedDraft, outX, outY, xCell, yCell)

outFile = netCDF4.Dataset(outFileName, 'w', format='NETCDF3_64BIT_OFFSET')
outFile.createDimension('Time', size=None)
outFile.createDimension('nCells', size=nCells)
outFile.createDimension('StrLen', size=StrLen)

outVar = outFile.createVariable('xtime', 'S1', ('Time', 'StrLen'))
for tIndex in range(2):
    outVar[tIndex, :] = netCDF4.stringtochar(xtime[tIndex])

outVar.setncatts({'units': 'unitless'})
writeVar('landIcePressureForcing', landIcePressure, 
         {'units': 'Pa', 
          'long_name': 'Pressure defined at the sea surface due to land ice'})
writeVar('landIceFractionForcing', landIceFraction, 
         {'units': 'unitless', 
          'long_name': 'The fraction of each cell covered by land ice'})
writeVar('landIceDraftForcing', landIceDraft, 
         {'units': 'unitless', 
          'long_name': 'The elevation of the interface between land ice and '
                       'the ocean'})

outFile.close()

