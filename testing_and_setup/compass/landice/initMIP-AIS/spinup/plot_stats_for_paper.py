#!/usr/bin/env python
'''
Script to plot time-series from globalStats file for AIS 20km spinup.
For MALI description paper.
'''

import sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

rhoi = 910.0

print "** Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser(description=__doc__)
parser.add_option("-1", dest="file1inName", help="input filename", default="globalStats.nc", metavar="FILENAME")
parser.add_option("-2", dest="file2inName", help="input filename", metavar="FILENAME")
parser.add_option("-3", dest="file3inName", help="input filename", metavar="FILENAME")
options, args = parser.parse_args()


# create axes to plot into
fig = plt.figure(1, facecolor='w')

nrow=2
ncol=2


axVAF = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('VAF (Gt)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()
axX = axVAF

axGrdArea = fig.add_subplot(nrow, ncol, 2, sharex=axX)
plt.xlabel('Year')
plt.ylabel('grounded area (km$^2$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axFltArea = fig.add_subplot(nrow, ncol, 4, sharex=axX)
plt.xlabel('Year')
plt.ylabel('floating area (km$^2$)')
axFltArea.yaxis.get_major_formatter().set_powerlimits((0,1))
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axGLflux = fig.add_subplot(nrow, ncol, 3, sharex=axX)
plt.xlabel('Year')
plt.ylabel('GL flux (Gt yr$^{-1}$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

#axCalvFlux = fig.add_subplot(nrow, ncol, 8, sharex=axX)
#plt.xlabel('Year')
#plt.ylabel('calving flux (kg/yr)')
##plt.xticks(np.arange(22)*xtickSpacing)
#plt.grid()



def plotStat(fname):
    print "Reading and plotting file: " + fname

    name = fname

    s=1
    e=200+1
    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][s:e]/365.0

#    vol = f.variables['totalIceVolume'][:]
#    axVol.plot(yr, vol, label=name)

    VAF = f.variables['volumeAboveFloatation'][s:e] / 1.0e12 * rhoi
    axVAF.plot(yr, VAF, label=name)

#    volGround = f.variables['groundedIceVolume'][:]
#    axVolGround.plot(yr, volGround, label=name)

#    volFloat = f.variables['floatingIceVolume'][:]
#    axVolFloat.plot(yr, volFloat, label=name)

    areaGrd = f.variables['groundedIceArea'][s:e] / 1000.0 / 1000.0
    axGrdArea.plot(yr, areaGrd, label=name)

    areaFlt = f.variables['floatingIceArea'][s:e] / 1000.0 / 1000.0
    axFltArea.plot(yr, areaFlt, label=name)

    GLflux = f.variables['groundingLineFlux'][s:e] / 1.0e12
    axGLflux.plot(yr, GLflux, label=name)

#    calvFlux = f.variables['totalCalvingFlux'][:]
#    axCalvFlux.plot(yr, calvFlux, label=name)


    f.close()


plotStat(options.file1inName)

if(options.file2inName):
   plotStat(options.file2inName)

if(options.file3inName):
   plotStat(options.file3inName)




print "Generating plot."

fig.tight_layout()
plt.show()
