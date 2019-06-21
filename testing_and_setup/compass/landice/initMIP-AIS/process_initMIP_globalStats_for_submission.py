#!/usr/bin/env python
'''
This script takes a globalStats.nc output file and processes it to follow the initMIP-AIS 
requirements for submission of scalar variables described here:
        http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica
Note that because fluxes at time 0 are meant to be from the end of the spinup,
a globalStats file from the final year of the spinup is also required.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser


print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n")
parser = OptionParser()
parser.description = __doc__
parser.add_option("-f", "--file", dest="inputFile", help="name of source (input) file.", default="globalStats.nc", metavar="FILENAME")
parser.add_option("-s", "--spinup", dest="spinupFile", help="name of file with spinup data.  Should only include the final year of spinup!  Values are evenly averaged.", default="spinupFinalYear.nc", metavar="FILENAME")
parser.add_option("-o", "--output", dest="outputFile", help="name of output file.", default="initMIP_scalar_stats.nc", metavar="FILENAME")
for option in parser.option_list:
   if option.default != ("NO", "DEFAULT"):
      option.help += (" " if option.help else "") + "[default: %default]"
      options, args = parser.parse_args()


rhoi = 910.0
print("Using ice density of:".format(rhoi))


inputData = Dataset(options.inputFile,'r')
years = inputData.variables['daysSinceStart'][:]/365.0
startYr = np.floor(years).min()
endYr = np.floor(years).max()
yrList = np.arange(startYr, endYr+1)

dt = inputData.variables['deltat'][:]
secInYr = 3600.0*24.0*365.0

evenYrInd = np.zeros_like(yrList).astype(int)
for t in range(len(yrList)):
   yr = yrList[t]
   evenYrInd[t] = np.nonzero(years >= yr)[0][0]  # get first index after the even year (just in case we don't exactly have an even year)
   print("Found time "+ str(years[evenYrInd[t]]) + " for index "+ str(t))

# also need spinup for flux data at t=0 :(
spinupData = Dataset(options.spinupFile,'r')




def annualAverage(varName):
   newVar = np.zeros_like(yrList)

   varSpinupData = spinupData.variables[varName][:]
   varData = inputData.variables[varName][:]
   for t in range(len(yrList)):
      if t == 0:
         newVar[t] = varSpinupData.mean()
      else:
         yr = yrList[t]
         #print yr
         ind = np.nonzero(np.logical_and(years>yr-1, years<=yr))[0]  # todo could skip redoing this for each var if slow
         assert dt[ind].sum() == secInYr, "sec in yr={}".format(dt[ind].sum())
         newVar[t] = (varData[ind] * dt[ind]).sum() / secInYr
   #print newVar
   return newVar

#varsToSkip = ('xtime','daysSinceStart','deltat')
#for var in inputData.variables:
#   if var in varsToSkip:
#      continue  # skip the rest of this iteration
#   print "Processing :", var
#   varData = inputData.variables[var][:]
#   newVar = np.zeros_like(yrList)
#   for t in range(len(yrList)):
#      yr = yrList[t]
#      ind = np.nonzero(np.logical_and(years>yr-1, years<=yr))[0]  # todo could skip redoing this for each var if slow
#      assert dt[ind].sum() == secInYr
#      newVar[t] = (varData[ind] * dt[ind]).sum() / secInYr
#   print newVar




# Create file and build time dimension
fout = Dataset(options.outputFile, 'w')

fout.createDimension('time', len(yrList))
timeVar = fout.createVariable('time', 'f', ('time',))
timeVar.units = "seconds since 2007-01-01 00:00:00"
timeVar.standard_name = "model_time"
timeVar[:] = yrList * secInYr

# process variables required
# Note: instructions request single precision

# land ice mass
lim = inputData.variables['totalIceVolume'][:] * rhoi
limVar = fout.createVariable('lim','f', ('time',))
limVar.units = "kg"
limVar.standard_name = "land_ice_mass"
limVar[:] = lim[evenYrInd]

# VAF
limnsw = inputData.variables['volumeAboveFloatation'][:] * rhoi
limnswVar = fout.createVariable('limnsw','f', ('time',))
limnswVar.units = "kg"
limnswVar.standard_name = "land_ice_mass_not_displacing_sea_water"
limnswVar[:] = limnsw[evenYrInd]

# grounded area
iareag = inputData.variables['groundedIceArea'][:]
var = fout.createVariable('iareag','f', ('time',))
var.units = "m^2"
var.standard_name = "grounded_land_ice_area"
var[:] = iareag[evenYrInd]

# floating area
iareaf = inputData.variables['floatingIceArea'][:]
var = fout.createVariable('iareaf','f', ('time',))
var.units = "m^2"
var.standard_name = "floating_ice_shelf_area"
var[:] = iareaf[evenYrInd]

# -- flux fields --

# SMB
var = fout.createVariable('tendacabf','f', ('time',))
var.units = "kg s^{-1}"
var.standard_name = "tendency_of_land_ice_mass_due_to_surface_mass_balance"
var[:] = annualAverage('totalSfcMassBal') / secInYr

# BMB
var = fout.createVariable('tendlibmassbf','f', ('time',))
var.units = "kg s^{-1}"
var.standard_name = "tendency_of_land_ice_mass_due_to_basal_mass_balance"
var[:] = annualAverage('totalBasalMassBal') / secInYr

# calving
var = fout.createVariable('tendlicalvf','f', ('time',))
var.units = "kg s^{-1}"
var.standard_name = "tendency_of_land_ice_mass_due_to_calving"
var[:] = -1.0 * annualAverage('totalCalvingFlux') / secInYr

# GL flux
var = fout.createVariable('tendligroundf','f', ('time',))
var.units = "kg s^{-1}"
var.standard_name = "tendency_of_grounded_ice_mass"
var[:] = -1.0 * annualAverage('groundingLineFlux') / secInYr

fout.Author="Matthew Hoffman (mhoffman@lanl.gov)"
fout.Model="MALI (MPAS-Albany Land Ice)"
fout.Variables="Scalar variables"
fout.Notes="Experiments performed at Los Alamos National Laboratory using the Edison supercomputer at National Energy Research Scientific Computing Center at Lawrence Berkeley National Laboratory.  Experiments performed by Matthew Hoffman, Tong Zhang, Stephen Price, and Mauro Perego."
fout.Date="28-Aug-2018"


fout.close()
inputData.close()
spinupData.close()


print("Complete.")

