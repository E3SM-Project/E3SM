#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
from netCDF4 import Dataset
import subprocess
import glob
import os

in_forcing_file = "forcing_data_init.nc"
out_forcing_file = "forcing_data_updated.nc"
out_forcing_link = "forcing_data.nc"

subprocess.check_call(['cp', in_forcing_file, out_forcing_file])
subprocess.check_call(['ln', '-sfn', out_forcing_file, out_forcing_link])

lastFileName = sorted(glob.glob('timeSeriesStatsMonthly*.nc'))[-1]

outFile = Dataset(out_forcing_file, 'r+')

inFile = Dataset('init.nc', 'r')
areaCell = inFile.variables['areaCell'][:]
ssh0 = inFile.variables['ssh'][0, :]
nCells = len(inFile.dimensions['nCells'])
inFile.close()

evaporationFluxVar = outFile.variables['evaporationFlux']
seaIceSalinityFluxVar = outFile.variables['seaIceSalinityFlux']
seaIceHeatFluxVar = outFile.variables['seaIceHeatFlux']
evaporationFlux = evaporationFluxVar[0, :]

evapMask = evaporationFlux != 0.

evapArea = numpy.sum(areaCell*evapMask)

totalArea = numpy.sum(areaCell)

rho_sw = 1026.
cp_sw = 3.996e3
secPerYear = 365*24*60*60

sflux_factor = 1.0
hflux_factor = 1.0/(rho_sw*cp_sw)

Tsurf = -1.9
Ssurf = 33.8

meanSSH0 = numpy.sum(ssh0*evapMask*areaCell)/evapArea

inFile = Dataset(lastFileName, 'r')
ssh = inFile.variables['timeMonthly_avg_ssh'][0, :]
inFile.close()
meanSSH = numpy.sum(ssh*evapMask*areaCell)/evapArea

deltaSSH = max(meanSSH - meanSSH0, 0.)

# parameterize the outgoing flux as a "spillway" with the given width,
# chosen to be 500 m so that a 20-m excess height is dissipated in about 2 months
spillwayWidth = 500.
g = 9.81

flowRate = -spillwayWidth*numpy.sqrt(0.5*g*deltaSSH**3)

estimatedHeightChange = flowRate*30*24*60*60/totalArea

# evap (m/s) is only over evapArea and negative for evaporation rather than
# precipitation
meanEvapRate = flowRate/evapArea
evapRate = meanEvapRate*evapMask

evaporationFlux = evapRate*rho_sw

print('update evap: mean sea-level increase: {} m'.format(deltaSSH))

print('update evap: est. one-month reduction by evap: {} m'.format(
    estimatedHeightChange))

print('update evap: evaporation rate: {} m/yr'.format(meanEvapRate*secPerYear))

evaporationFluxVar[0, :] = evaporationFlux

seaIceSalinityFluxVar[0, :] = evapRate*Ssurf/sflux_factor

seaIceHeatFluxVar[0, :] = evapRate*Tsurf/hflux_factor

outFile.close()
