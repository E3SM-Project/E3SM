#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots Ice1r/Ice1ra for a series of resolutions.

See bisicles_ssa_tsai_datasheet.pdf from Supplement from Asay-Davis paper for example.

by Matt Hoffman, modified from mismip+PlotGL.py script.
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

model = '_MPASLI'

def vafplot(resolution, color, marker):
    """
    add a plot of change in volume above floatation against time to current axes
    """

    fname = '{}/Ice1r/Ice1r{}.nc'.format(resolution, model)
    label = "{}".format(resolution)
    try:
       ncid = Dataset(fname, 'r')
       vaf1r = ncid.variables["iceVAF"][:]*(1e-3**3*1e-3)
       time1r = ncid.variables["time"][:]
       ncid.close()
    except:
       print "Failed to open file: {}. Skipping.".format(fname)
       vaf1r = []
       time1r = []
       
    # repeat for 1ra
    fname = '{}/Ice1ra/Ice1ra{}.nc'.format(resolution, model)
    label = "{}".format(resolution)
    try:
       ncid = Dataset(fname, 'r')
       vaf1ra = ncid.variables["iceVAF"][:]*(1e-3**3*1e-3)
       time1ra = ncid.variables["time"][:]
       ncid.close()
    except:
       print "Failed to open file: {}. Skipping.".format(fname)
       vaf1ra = []
       time1ra = []

    # repeat for 1rax
    fname = '{}/Ice1rax/Ice1rax{}.nc'.format(resolution, model)
    label = "{}".format(resolution)
    try:
       ncid = Dataset(fname, 'r')
       vaf1rax = ncid.variables["iceVAF"][:]*(1e-3**3*1e-3)
       time1rax = ncid.variables["time"][:]
       ncid.close()
    except:
       print "Failed to open file: {}. Skipping.".format(fname)
       vaf1rax = []
       time1rax = []


    vaf = np.concatenate((vaf1r, vaf1ra, vaf1rax))
    time = np.concatenate((time1r, time1ra, time1rax))
    vafDelta = vaf - vaf[0]

    plt.plot(time, vafDelta, 'o-', mfc=color, ms=4,
                 color='black', label=label, marker=marker)
    return np.max(vafDelta)

def glplot(ncfile, times, colora, label):
    """
    add a plot of grounding line points to current axes.
    makes use of the numpy.ma.MaskedArray when reading xGL,yGL
    """
    try:
       ncid = Dataset(ncfile, 'r')
    except:
       print "Failed to open file: {}. Skipping.".format(ncfile)
       return 350.0, 500.0

    time = ncid.variables["time"][:]
    lxmax = 0.0
    lxmin = 800.0
    for i in range(0, len(times)):
        seq = (time == times[i])
        xGL = ncid.variables["xGL"][:, seq]*1e-3
        lxmax = max(np.max(xGL), lxmax)
        lxmin = min(np.min(xGL), lxmin)
        yGL = ncid.variables["yGL"][:, seq]*1e-3
        plt.plot(xGL, yGL, 's', ms=2, mfc=colora[i],
                 mec=colora[i], label=label + ', t = ' + format(times[i]))
    return lxmin, lxmax

plt.figure(figsize=(12, 5))

plt.subplot(121)

times=[0, 100]
xmin, xmax = glplot('4000m/Ice1r/Ice1r' + model + '.nc', times, ['black', 'black'], '4000m')
xmin, xmax = glplot('2000m/Ice1r/Ice1r' + model + '.nc', times, ['blue', 'blue'], '2000m')
xmin, xmax = glplot('1000m/Ice1r/Ice1r' + model + '.nc', times, ['purple', 'purple'], '1000m')
xmin, xmax = glplot('500m/Ice1r/Ice1r' + model + '.nc', times, ['red', 'red'], '500m')
xmin, xmax = glplot('250m/Ice1r/Ice1r' + model + '.nc', times, ['orange', 'orange'], '250m')
#plt.legend(frameon=True, borderaxespad=0, loc='right')
plt.xlabel(r'$x$ (km)')
plt.ylabel(r'$y$ (km)')


# ======================================
#ax = plt.figure(figsize=(7, 5))
plt.subplot(122)
yrange=[-5, 1]
plt.plot([100, 100], yrange, color="grey")
plt.plot([200, 200], yrange, color="grey")
plt.xlim([0, 500])
plt.ylim(yrange)

xtlocs = [0, 100, 200, 300, 400, 500]
plt.xticks(xtlocs, xtlocs)
plt.xlabel(r'Time,  $t$ (a)')
plt.ylabel(r'$\Delta$ VAF (1000 km$^3$)')

maxa = vafplot('4000m', 'black', 'o')

maxa = vafplot('2000m', 'blue', 'o')

maxa = vafplot('1000m', 'purple', 'o')

maxa = vafplot('500m', 'red', 'o')

maxa = vafplot('250m', 'orange', 'o')

plt.legend(loc='upper right', ncol=1, frameon=True, borderaxespad=0)

plotname = 'mismip+ResolutionAnalysis.pdf'
plt.savefig(plotname)
plt.show()
plt.close()

print 'Created test plot', plotname
