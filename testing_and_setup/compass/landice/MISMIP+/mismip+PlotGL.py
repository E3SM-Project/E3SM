#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot grounded area against time for all 3 MISMIP+ experiments (7 branches)
plot the grounding lines at 0,  100, 200 years

Created on Tue Sep  8 09:29:09 2015

@author: s.l.cornford@bris.ac.uk

Modified by William Lipscomb
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

model = '_MPASLI'

def tscale(time):
    """
    scale time to sqrt(time) to emphasize earlier times
    """
    return np.sqrt(time)

def intscale(time):
    """
    inverse of tscale
    """
    return time**2

def garplot(ncfile, label, color, marker):
    """
    add a plot of grounded area against time to current axes
    """
    ncid = Dataset(ncfile, 'r')
    gar = ncid.variables["groundedArea"][:]*1e-6*1e-3
    time = ncid.variables["time"][:]
    if label == 'nolabel':
        plt.plot(tscale(time), gar, 'o-', mfc=color,
                 color='black', label="", marker=marker)
    else:    
        plt.plot(tscale(time), gar, 'o-', mfc=color,
                 color='black', label=label, marker=marker)
    ncid.close()
    return np.max(gar)

def glplot(ncfile, times, colora, label):
    """
    add a plot of grounding line points to current axes.
    makes use of the numpy.ma.MaskedArray when reading xGL,yGL
    """
    ncid = Dataset(ncfile, 'r')
    time = ncid.variables["time"][:]
    lxmax = 0.0
    lxmin = 800.0
    for i in range(0, len(times)):
        seq = (time == times[i])
        xGL = ncid.variables["xGL"][:, seq]*1e-3
        lxmax = max(np.max(xGL), lxmax)
        lxmin = min(np.min(xGL), lxmin)
        yGL = ncid.variables["yGL"][:, seq]*1e-3
        plt.plot(xGL, yGL, 's', ms=3, mfc=colora[i],
                 mec=colora[i], label=label + ', t = ' + format(times[i]))
    return lxmin, lxmax

plt.figure(figsize=(7, 10))

plt.subplot(211)

xmin, xmax = glplot('Ice1r/Ice1r' + model + '.nc', [0, 100], ['black', 'red'], 'Ice1r')
xmaxplot = xmax  # xmax based on Ice1r at t = 0
xmin, xmax = glplot('Ice1rr/Ice1rr' + model + '.nc', [200], ['yellow'], 'Ice1rr')
xminplot = xmin  # xmin based on Ice1rr at t = 200
plt.xlim([xminplot-20.0, xmaxplot+50.0])
xmin, xmax = glplot('Ice1ra/Ice1ra' + model + '.nc', [200], ['orange'], 'Ice1ra')
xmin, xmax = glplot('Ice2r/Ice2r' + model + '.nc', [100], ['blue'], 'Ice2r')
xmin, xmax = glplot('Ice2rr/Ice2rr' + model + '.nc', [200], ['pink'], 'Ice2rr')
xmin, xmax = glplot('Ice2ra/Ice2ra' + model + '.nc', [200], ['purple'], 'Ice2ra')

plt.legend(frameon=True, borderaxespad=0, loc='right')
plt.xlabel(r'$x$ (km)')
plt.ylabel(r'$y$ (km)')

#ax = plt.figure(figsize=(7, 5))
plt.subplot(212)
plt.plot(tscale([100, 100]), [0, 100], color="grey")
plt.plot(tscale([200, 200]), [0, 100], color="grey")
plt.xlim(tscale([0, 1000]))
plt.ylim([25, 40])

xtlocs = tscale([0, 10, 50, 100, 200, 400, 800])
plt.xticks(xtlocs, intscale(xtlocs))
plt.xlabel(r'Time,  $t$ (a)')
plt.ylabel(r'Grounded area (1000 km$^3$)')

#Ice0
maxa = garplot('Ice0/Ice0' + model + '.nc', 'Ice0', 'grey', 'd')

#Ice1
maxa = garplot('Ice1r/Ice1r' + model + '.nc', 'Ice1r', 'red', 'o')
maxa = garplot('Ice1rr/Ice1rr' + model + '.nc', 'Ice1rr', 'purple', 'o')
maxa = garplot('Ice1ra/Ice1ra' + model + '.nc', 'Ice1ra', 'orange', 'o')
maxa = garplot('Ice1rrx/Ice1rrx' + model + '.nc', 'nolabel', 'purple', 'o')
maxa = garplot('Ice1rax/Ice1rax' + model + '.nc', 'nolabel', 'orange', 'o')

#Ice2
maxa = garplot('Ice2r/Ice2r' + model + '.nc', 'Ice2r', 'blue', 's')
maxa = garplot('Ice2rr/Ice2rr' + model + '.nc', 'Ice2rr', 'pink', 's')
maxa = garplot('Ice2ra/Ice2ra' + model + '.nc', 'Ice2ra', 'yellow', 's')
maxa = garplot('Ice2rrx/Ice2rrx' + model + '.nc', 'nolabel', 'pink', 's')
maxa = garplot('Ice2rax/Ice2rax' + model + '.nc', 'nolabel', 'yellow', 's')

plt.legend(loc='lower left', ncol=2, frameon=True, borderaxespad=0)

plotname = 'mismip+PlotGL.pdf'
plt.savefig(plotname)
plt.close()

print 'Created test plot', plotname
