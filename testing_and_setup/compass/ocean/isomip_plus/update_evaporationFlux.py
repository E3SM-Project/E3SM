#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
from netCDF4 import Dataset

from optparse import OptionParser

import subprocess

import glob

parser = OptionParser()

parser.add_option("--in_forcing_file", type="string",
                  default="forcing_data_init.nc", dest="in_forcing_file")
parser.add_option("--out_forcing_file", type="string",
                  default="forcing_data_updated.nc", dest="out_forcing_file")
parser.add_option("--out_forcing_link", type="string",
                  default="forcing_data.nc", dest="out_forcing_link")
parser.add_option("--avg_months", type="int", default=3, dest="avg_months")

options, args = parser.parse_args()

subprocess.check_call(['cp',
                       options.in_forcing_file, options.out_forcing_file])
subprocess.check_call(['ln', '-sfn',
                       options.out_forcing_file, options.out_forcing_link])

fluxFiles = sorted(glob.glob('timeSeriesStatsMonthly*.nc'))
index = max(0, len(fluxFiles)-options.avg_months)
print(index)
fluxFiles = fluxFiles[index:]

outFile = Dataset(options.out_forcing_file, 'r+')

inFile = Dataset('init.nc', 'r')
areaCell = inFile.variables['areaCell'][:]
fraction = inFile.variables['landIceFraction'][0, :]
meanIceArea = numpy.sum(fraction*areaCell)

nCells = len(inFile.dimensions['nCells'])
inFile.close()

rho_sw = 1026.
cp_sw = 3.996e3
secPerYear = 365*24*60*60

sflux_factor = 1.0
hflux_factor = 1.0/(rho_sw*cp_sw)

Tsurf = -1.9
Ssurf = 33.8

meanMeltFlux = 0.0
for fileName in fluxFiles:
    inFile = Dataset(fileName, 'r')
    freshwaterFlux = \
        inFile.variables['timeMonthly_avg_landIceFreshwaterFlux'][0, :]
    meanMeltFlux += numpy.sum(freshwaterFlux*areaCell)
    inFile.close()

meanMeltFlux /= len(fluxFiles)

# convert to volume flux in m^3/s
meanMeltFlux /= rho_sw

meanMeltRate = meanMeltFlux*secPerYear/meanIceArea
print('mean melt rate: {} m/yr'.format(meanMeltRate))

area = numpy.sum(areaCell)

meanSeaLevelRiseRate = meanMeltFlux*secPerYear/area
print('mean rate of sea-level change: {} m/yr'.format(meanSeaLevelRiseRate))

evaporationFluxVar = outFile.variables['evaporationFlux']
seaIceSalinityFluxVar = outFile.variables['seaIceSalinityFlux']
seaIceHeatFluxVar = outFile.variables['seaIceHeatFlux']
evaporationFlux = evaporationFluxVar[0, :]

evapMask = evaporationFlux != 0.

evapArea = numpy.sum(areaCell*evapMask)

evapRate = -meanMeltFlux/evapArea

print('evaporation rate: {} m/yr'.format(evapRate*secPerYear))

evaporationFlux[evapMask] = evapRate*rho_sw
evaporationFluxVar[0, :] = evaporationFlux

flux = seaIceSalinityFluxVar[0, :]
flux[evapMask] = evapRate*Ssurf/sflux_factor
seaIceSalinityFluxVar[0, :] = flux

flux = seaIceHeatFluxVar[0, :]
flux[evapMask] = evapRate*Tsurf/hflux_factor
seaIceHeatFluxVar[0, :] = flux

outFile.close()
