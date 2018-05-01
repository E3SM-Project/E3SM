#!/usr/bin/env python
import numpy
from netCDF4 import Dataset

from optparse import OptionParser

import os
import os.path

import subprocess

parser = OptionParser()

parser.add_option("--in_fluxes_file", type="string", default="land_ice_fluxes.nc", dest="in_fluxes_file")
parser.add_option("--in_forcing_file", type="string", default="forcing_data_init.nc", dest="in_forcing_file")
parser.add_option("--out_forcing_file", type="string", default="forcing_data_updated.nc", dest="out_forcing_file")
parser.add_option("--out_forcing_link", type="string", default="forcing_data.nc", dest="out_forcing_link")
parser.add_option("--avg_years", type="float", default=0.25, dest="avg_years")

options, args = parser.parse_args()

subprocess.check_call(['cp', options.in_forcing_file, options.out_forcing_file])
subprocess.check_call(['ln', '-sfn', options.out_forcing_file, options.out_forcing_link])

inFile = Dataset(options.in_fluxes_file,'r')
outFile = Dataset(options.out_forcing_file,'r+')

areaCell = inFile.variables['areaCell'][:]
times = inFile.variables['daysSinceStartOfSim'][:]/365.

nCells = len(inFile.dimensions['nCells'])
nTime = len(inFile.dimensions['Time'])

tIndices = numpy.nonzero(times >= times[-1]-options.avg_years)[0]

print len(tIndices)

rho_sw = 1026.
cp_sw = 3.996e3
secPerYear = 365*24*60*60

sflux_factor = 1.0
hflux_factor = 1.0/(rho_sw*cp_sw)

Tsurf = -1.9
Ssurf = 33.8

meanMeltFlux = 0.0
meanIceArea = 0.0
for tIndex in tIndices:
  freshwaterFlux = inFile.variables['landIceFreshwaterFlux'][tIndex,:]
  fraction = inFile.variables['landIceFraction'][tIndex,:]
  meanMeltFlux += numpy.sum(freshwaterFlux*areaCell)
  meanIceArea += numpy.sum(fraction*areaCell)

meanMeltFlux /= len(tIndices)
meanIceArea /= len(tIndices)

# convert to volume flux in m^3/s
meanMeltFlux /= rho_sw

meanMeltRate = meanMeltFlux*secPerYear/meanIceArea
print 'mean melt rate:', meanMeltRate, 'm/yr'

area = numpy.sum(areaCell)

meanSeaLevelRiseRate = meanMeltFlux*secPerYear/area
print 'mean rate of sea-level change:', meanSeaLevelRiseRate, 'm/yr'

evaporationFluxVar = outFile.variables['evaporationFlux']
seaIceSalinityFluxVar = outFile.variables['seaIceSalinityFlux']
seaIceHeatFluxVar = outFile.variables['seaIceHeatFlux']
evaporationFlux = evaporationFluxVar[0,:]

evapMask = evaporationFlux != 0.

evapArea = numpy.sum(areaCell*evapMask)

evapRate = -meanMeltFlux/evapArea

print 'evaporation rate:', evapRate*secPerYear, 'm/yr'

evaporationFlux[evapMask] = evapRate*rho_sw
evaporationFluxVar[0,:] = evaporationFlux

flux = seaIceSalinityFluxVar[0,:]
flux[evapMask] = evapRate*Ssurf/sflux_factor
seaIceSalinityFluxVar[0,:] = flux

flux = seaIceHeatFluxVar[0,:]
flux[evapMask] = evapRate*Tsurf/hflux_factor
seaIceHeatFluxVar[0,:] = flux

inFile.close()
outFile.close()
