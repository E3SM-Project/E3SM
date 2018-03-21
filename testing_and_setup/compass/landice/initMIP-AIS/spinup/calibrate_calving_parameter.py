#!/usr/bin/env python
'''
Script to calibrate calving parameter for eigencalving from an optimized model state.
'''

import sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
from scipy import interpolate



print "** Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.add_option("-i", "--in", dest="fileinName", help="input filename.", default="landice_grid.nc", metavar="FILENAME")
parser.add_option("-o", "--out", dest="fileoutName", help="output filename.", default="eigencalvingParameter.nc", metavar="FILENAME")
options, args = parser.parse_args()

time=-1 # use last time in file

# ----
# Get the input file
filein = Dataset(options.fileinName,'r')
# get needed vars from input
xCell = filein.variables['xCell'][:]
yCell = filein.variables['yCell'][:]
cellsOnCell = filein.variables['cellsOnCell'][:]
nEdgesOnCell = filein.variables['nEdgesOnCell'][:]
eMax = filein.variables['eMax'][time,:]
eMin = filein.variables['eMin'][time,:]
speed = filein.variables['surfaceSpeed'][time,:]
cellMask = filein.variables['cellMask'][time,:]
thickness = filein.variables['thickness'][time,:]
bedTopography = filein.variables['bedTopography'][time,:]
nCells = len(filein.dimensions['nCells'])



# ----
# Define the new file to be output - will be clobbered!
fileout = Dataset(options.fileoutName,"w",format=filein.file_format)
# create field to be populated
datatype = filein.variables['thickness'].dtype  # Get the datatype for float
fileout.createDimension('Time', 1)
fileout.createDimension('nCells', nCells)
fileout.createVariable('eigencalvingParameter', datatype, ('Time', 'nCells',))

# calc K
eps = 1.0e-30

# make mask
calvMask = np.logical_or(eMax<=0.0, eMin<=0.0)  # mask out compressional strain
calvMask = np.logical_or(calvMask, np.invert((cellMask&4).astype(bool)) )   # mask out grounded ice  (not floating)
calvMask = np.logical_or(calvMask, np.invert((cellMask&16).astype(bool)) )   # mask out ice not at the dynamic margin
#K = calvMask  # use this to print out mask

# calculate K
#K = speed / (np.maximum(0.0, eMax) * np.maximum(0.0, eMin) + eps)  # version with MA
K = ma.array(speed / (np.maximum(0.0, eMax) * np.maximum(0.0, eMin)), mask=calvMask)

K2 = K.filled(0)  # get version with 0s in masked area
K2 = np.clip(K2, 1.0e15, 1.0e17) # clip extremes from range

calvMask2 = np.logical_or(calvMask, K2>1.0e17)
calvMask2 = np.logical_or(calvMask2, K2<1.0e15)
K3 = ma.array(K2, mask=calvMask2)

#K2=K2*1.5  # increase it, seems to small

method = 2

if method == 1:
   # calibrate locally then extrapolate
   print "Start extrapolation!"
   filledCells = np.invert(calvMask)  # init this mask to where we have valid K
   print "nCells=", nCells
   print "{} cells left for extrapolation.".format(nCells - np.count_nonzero(filledCells))
   while np.count_nonzero(filledCells) != nCells:
      filledCellsNew = np.copy(filledCells)
      searchCells = np.where(filledCells==0)[0]
      #print "searching:", searchCells
      for iCell in searchCells:
         neighbors = cellsOnCell[iCell, :nEdgesOnCell[iCell]]-1
         goodNeigh = neighbors[np.nonzero(filledCells[neighbors])]
         if len(goodNeigh>0):
           K2[iCell] = K2[goodNeigh].mean()
           #K2[iCell] = 10.0**(np.log10(K2[goodNeigh]).mean())
           filledCellsNew[iCell] = True
      filledCells = np.copy(filledCellsNew)  # update mask
      print "{} cells left for extrapolation.".format(nCells - np.count_nonzero(filledCells))


elif method == 2:
   # set by region
   fregion =  Dataset('regionMasks.nc','r')
   nRegions = len(fregion.dimensions['nRegions'])
   names = fregion.variables['regionNames'][:]
   regions = fregion.variables['regionCellMasks'][:]
   megaRegions={   'RF':     {'imbie':[1,2,3,27], 'k':3.7e16},
                   'penin':  {'imbie':[23, 24, 25, 26], 'k':2.0e16},
                   'amund': {'imbie':[22, 21], 'k':2.0e16},
                   'ross': {'imbie':[20,19,18,17,16], 'k':3.5e16},
                   'eais': {'imbie':[15,14,13,12,7,6,5,4], 'k':0.3e16},
                   'amery': {'imbie':[8,9,10,11], 'k':0.3e16}
               }
   imbieNums=np.zeros((nRegions,))
   megaregionMean = np.zeros((len(megaRegions),))
   for r in range(nRegions):
      imbieNums[r] = int("".join(names[r, 16:18]))
   i=0
   for m in megaRegions:
      allInd=np.zeros((1,), dtype=int)
      for r in range(nRegions):
         if imbieNums[r] in megaRegions[m]['imbie']:
            thisRegionInd = np.nonzero(regions[:,r] == 1)[0]
            allInd = np.append(allInd, thisRegionInd)

      allInd2 = np.sort(allInd[1:])
      print len(allInd2)
      megaregionMean[i] = K3[allInd2].mean()
      K2[allInd2] = megaregionMean[i]
      K2[allInd2] = megaRegions[m]['k']
      i = i + 1

   # there is some places in no region :(
   ind = np.nonzero(regions.sum(axis=1)==0)[0]
   K2[ind] = K2[ np.nonzero(regions[:,4]==1)[0][0]  ]   # sticking in value from imbie 4

#fileout.variables['eigencalvingParameter'][0,:] = K2
##fileout.variables['eigencalvingParameter'][0,:] = K.filled(0)
#fileout.close()
#sys.exit()




print ""

print "Adding buffer increase"
filledCells = (thickness>0.0)
print "{} cells left for buffer increase.".format(nCells - np.count_nonzero(filledCells))
while np.count_nonzero(filledCells) != nCells:
   filledCellsNew = np.copy(filledCells)
   searchCells = np.where(filledCells==0)[0]
   #print "searching:", searchCells
   for iCell in searchCells:
      neighbors = cellsOnCell[iCell, :nEdgesOnCell[iCell]]-1
      goodNeigh = neighbors[np.nonzero(filledCells[neighbors])]
      if len(goodNeigh>0):
        K2[iCell] = K2[goodNeigh].mean() * 1.01   # increase K in each layer beyond the initial ice by this factor
        filledCellsNew[iCell] = True
   filledCells = np.copy(filledCellsNew)  # update mask
   print "{} cells left for buffer increase.".format(nCells - np.count_nonzero(filledCells))


# add huge values at edge of mesh (two rows)
filledCells = filledCells*0
for iCell in range(nCells):
   if cellsOnCell[iCell, :nEdgesOnCell[iCell]].min() == 0:
      filledCells[iCell] = 1
      K2[iCell] = 1.0e17
searchCells = np.where(filledCells==0)[0]
for iCell in searchCells:
      neighbors = cellsOnCell[iCell, :nEdgesOnCell[iCell]]-1
      goodNeigh = neighbors[np.nonzero(filledCells[neighbors])]
      if len(goodNeigh>0):
         K2[iCell] = 1.0e17

#indGood = np.where(calvMask == 0)[0]
#myInterp = interpolate.LinearNDInterpolator(np.stack((xCell[indGood], yCell[indGood])).T, K2[indGood])
#indBad = np.where(calvMask == 1)[0]
#K2[indBad] = myInterp(np.stack((xCell[indBad], yCell[indBad])).T)

fileout.variables['eigencalvingParameter'][0,:] = K2
#fileout.variables['eigencalvingParameter'][0,:] = K.filled(0)
fileout.close()

print "Wrote field to " + options.fileoutName
