#!/usr/bin/env python
'''
Script to calculate annual time postings for selected regional stats output.
For flux fields calculates time averages, for state fields copies snapshots.
Output file is in MALI regionalStats format.
Note that because fluxes at time 0 are meant to be from the end of the spinup,
a globalStats file from the final year of the spinup is also required.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import sys


print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n")
parser = OptionParser()
parser.description = __doc__
parser.add_option("-f", "--file", dest="inputFile", help="name of source (input) file.", default="regionalStats.nc", metavar="FILENAME")
parser.add_option("-s", "--spinup", dest="spinupFile", help="name of file with spinup data.  Should only include the final year of spinup!  Values are evenly averaged.", default="spinupFinalYear.nc", metavar="FILENAME")
for option in parser.option_list:
   if option.default != ("NO", "DEFAULT"):
      option.help += (" " if option.help else "") + "[default: %default]"
      options, args = parser.parse_args()


rhoi = 910.0
print("Using ice density of: {}".format(rhoi))


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
   print("Found time " + str( years[evenYrInd[t]] ) + " for index " + str(t))

# also need spinup for flux data at t=0 :(
spinupData = Dataset(options.spinupFile,'r')




def annualAverage(varName):
   newVar = np.zeros((len(yrList), nRegions))

   varSpinupData = spinupData.variables[varName][:]
   varData = inputData.variables[varName][:]
   for t in range(len(yrList)):
      if t == 0:
         newVar[t,:] = varSpinupData.mean(axis=0)
      else:
         yr = yrList[t]
         #print yr
         ind = np.nonzero(np.logical_and(years>yr-1, years<=yr))[0]  # todo could skip redoing this for each var if slow
         assert dt[ind].sum() == secInYr, "sec in yr={}".format(dt[ind].sum())
         for r in range(nRegions):
            newVar[t,r] = (varData[ind,r] * dt[ind]).sum(axis=0) / secInYr
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
fout = Dataset('initMIP_scalar_regionalStats.nc','w')

fout.createDimension('Time', len(yrList))
timeVar = fout.createVariable('daysSinceStart', 'f', ('Time',))
timeVar[:] = yrList * 365.0

nRegions= len(inputData.dimensions['nRegions'])
fout.createDimension('nRegions', nRegions)

# process variables required
# Note: instructions request single precision

# land ice mass
lim = inputData.variables['regionalIceVolume'][:]
limVar = fout.createVariable('regionalIceVolume','f', ('Time','nRegions'))
limVar[:] = lim[evenYrInd,:]


# VAF
limnsw = inputData.variables['regionalVolumeAboveFloatation'][:]
limnswVar = fout.createVariable('regionalVolumeAboveFloatation','f', ('Time','nRegions'))
limnswVar[:] = limnsw[evenYrInd]

# grounded area
iareag = inputData.variables['regionalGroundedIceArea'][:]
var = fout.createVariable('regionalGroundedIceArea','f', ('Time','nRegions'))
var[:] = iareag[evenYrInd]

# floating area
iareaf = inputData.variables['regionalFloatingIceArea'][:]
var = fout.createVariable('regionalFloatingIceArea','f', ('Time','nRegions'))
var[:] = iareaf[evenYrInd]

# -- flux fields --


# BMB
var = fout.createVariable('regionalSumFloatingBasalMassBal','f', ('Time','nRegions'))
var[:] = annualAverage('regionalSumFloatingBasalMassBal')

# calving
var = fout.createVariable('regionalSumCalvingFlux','f', ('Time','nRegions'))
var[:] = annualAverage('regionalSumCalvingFlux')

# GL flux
var = fout.createVariable('regionalSumGroundingLineFlux','f', ('Time','nRegions'))
var[:] = annualAverage('regionalSumGroundingLineFlux')





fout.close()
inputData.close()
spinupData.close()


print("Complete.")

