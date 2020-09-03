#!/usr/bin/env python

from netCDF4 import Dataset as NetCDFFile

# Open the file, get needed dimensions
gridfile = NetCDFFile('./landice_grid.nc','r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
# Get variables
xCell = gridfile.variables['xCell']
yCell = gridfile.variables['yCell']
xEdge = gridfile.variables['xEdge']
yEdge = gridfile.variables['yEdge']
xVertex = gridfile.variables['xVertex']
yVertex = gridfile.variables['yVertex']
thickness = gridfile.variables['thickness']
bedTopography = gridfile.variables['bedTopography']
layerThicknessFractions = gridfile.variables['layerThicknessFractions']
SMB = gridfile.variables['sfcMassBal']
basalHeatFlux = gridfile.variables['basalHeatFlux']
surfaceAirTemperature = gridfile.variables['surfaceAirTemperature']
temperature = gridfile.variables['temperature']

thickness[:] = 200
bedTopography[:] = 0
basalHeatFlux[:] = 0.0
surfaceAirTemperature[:] = 270.15
temperature[:] = 270.15

gridfile.close()

