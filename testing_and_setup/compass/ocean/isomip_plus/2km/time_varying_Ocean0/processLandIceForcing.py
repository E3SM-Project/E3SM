#!/usr/bin/env python
import numpy
import netCDF4
import os
import sys


def writeVar(outVarName, field, attrs):
    outVar = outFile.createVariable(outVarName, 'f8', ('Time', 'nCells'))
    outVar[:, :] = field
    outVar.setncatts(attrs)


inMpasMeshFileName = 'init.nc'
outFileName = 'land_ice_forcing.nc'
initialScale = 0.1


if os.path.exists(outFileName):
    # nothing to do
    sys.exit(0)

# interpolate
mpasFile = netCDF4.Dataset(inMpasMeshFileName, 'r')
xCell = mpasFile.variables['xCell'][:]
yCell = mpasFile.variables['yCell'][:]
inLandIceDraft = mpasFile.variables['landIceDraft'][:]
inLandIcePressure = mpasFile.variables['landIcePressure'][:]
inLandIceFraction = mpasFile.variables['landIceFraction'][:]
mpasFile.close()

nCells = len(xCell)
StrLen = 64
nTime = 3
xtime = numpy.zeros((nTime), 'S64')

xtime[0] = "0001-01-01_00:00:00                                             "
xtime[1] = "0002-01-01_00:00:00                                             "
xtime[2] = "0003-01-01_00:00:00                                             "

landIcePressure = numpy.zeros((nTime, nCells), float)
landIcePressure[0, :] = inLandIcePressure
landIcePressure[1, :] = inLandIcePressure/initialScale
landIcePressure[2, :] = inLandIcePressure/initialScale
landIceFraction = numpy.zeros(landIcePressure.shape)
landIceFraction[0, :] = inLandIceFraction
landIceFraction[1, :] = inLandIceFraction
landIceFraction[2, :] = inLandIceFraction
landIceDraft = numpy.zeros(landIcePressure.shape)
landIceDraft[0, :] = inLandIceDraft
landIceDraft[1, :] = inLandIceDraft/initialScale
landIceDraft[2, :] = inLandIceDraft/initialScale

outFile = netCDF4.Dataset(outFileName, 'w', format='NETCDF3_64BIT_OFFSET')
outFile.createDimension('Time', size=None)
outFile.createDimension('nCells', size=nCells)
outFile.createDimension('StrLen', size=StrLen)

outVar = outFile.createVariable('xtime', 'S1', ('Time', 'StrLen'))
for tIndex in range(nTime):
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

