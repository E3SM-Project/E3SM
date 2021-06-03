#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
'''

import sys
import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

outfname = 'globalStats.nc'
runs=[ adir for adir in sorted(os.listdir('.')) if (os.path.isdir(adir) and 'adjust'in adir and os.path.isfile(os.path.join(adir, outfname)))]
print "Original run list:", runs

fig = plt.figure(1, facecolor='w')
nrow = 1; ncol = 1;
ax1 = fig.add_subplot(nrow, ncol, 1)
plt.ylabel('depth to mCDW (m)')
plt.xlabel('mean melt (m/yr)')
#plt.xticks(np.arange(30)*xtickSpacing)
plt.grid()

baseDepth = -700.0
runNumber = 0
for run in runs:
    print "Plotting run: " + run
    depthChange = float(run[6:])
    depth = baseDepth + depthChange

    f = netCDF4.Dataset(run+"/globalStats.nc", 'r')
    avgMelt = f.variables['avgSubshelfMelt'][-1]
    totalMelt = f.variables['totalFloatingBasalMassBal'][-1]

    ax1.plot(avgMelt, depth, '.k')
    print depthChange, depth, avgMelt

plt.show()
