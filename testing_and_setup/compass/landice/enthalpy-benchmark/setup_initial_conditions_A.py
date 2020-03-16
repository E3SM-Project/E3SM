#!/usr/bin/env python

from netCDF4 import Dataset as NetCDFFile

# Open the file, get needed dimensions
gridfile = NetCDFFile('./landice_grid.nc','r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
# Get variables
thickness = gridfile.variables['thickness']
bedTopography = gridfile.variables['bedTopography']
basalHeatFlux = gridfile.variables['basalHeatFlux']
surfaceAirTemperature = gridfile.variables['surfaceAirTemperature']
temperature = gridfile.variables['temperature']

thickness[:] = 1000
bedTopography[:] = 0
basalHeatFlux[:] = 0.042
surfaceAirTemperature[:] = 243.15
temperature[:] = 243.15

gridfile.close()

