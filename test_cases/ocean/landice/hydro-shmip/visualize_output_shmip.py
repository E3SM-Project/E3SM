#!/usr/bin/env python
'''
Plots profiles for hydro-ship test case
'''
import sys
import numpy as np
import netCDF4
#import datetime
# import math
# from pylab import *
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib import cm
# from matplotlib.contour import QuadContourSet
# import time

secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize", metavar="FILE")
parser.add_option("-t", "--time", dest="time", help="time step to visualize (0 based)", metavar="TIME")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")
parser.add_option("-3", dest="A3", action="store_true", help="plot GLADS results for experiment 3")
parser.add_option("-5", dest="A5", action="store_true", help="plot GLADS results for experiment 5")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

if not options.time:
	print "No time provided. Using time -1."
        time_slice = -1
else:
        time_slice = int(options.time)


if options.A3 and options.A5:
   sys.exit("Only one of -3 and -5 can be specified.")

f = netCDF4.Dataset(options.filename,'r')

#xtime = f.variables['xtime'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
yVertex = f.variables['yVertex'][:]
xEdge = f.variables['xEdge'][:]
#yEdge = f.variables['yEdge'][:]
h = f.variables['waterThickness'][time_slice,:]
u = f.variables['waterVelocityCellX'][time_slice,:]
N = f.variables['effectivePressure'][time_slice,:]
days = f.variables['daysSinceStart'][:]
#basalmelt = f.variables['basalMeltInput'][time_slice,:]
#surfmelt = f.variables['externalWaterInput'][time_slice,:]
xtime= f.variables['xtime']
areaCell = f.variables['areaCell'][:]

q = u*h


#print "attempting to get input data from landice_grid.nc!"
#fin = netCDF4.Dataset('landice_grid.nc','r')
#H = fin.variables['thickness'][0,:]

print "Using time level ", time_slice, ", which is xtime=",''.join( xtime[time_slice,:])

# Find center row  - currently files are set up to have central row at y=0
unique_ys=np.unique(yCell[:])
centerY=unique_ys[len(unique_ys)/2]
print "number of ys, center y index, center Y value", len(unique_ys), len(unique_ys)/2, centerY
ind = np.nonzero(yCell[:] == centerY)
x = xCell[ind]/1000.0

# calculate mean,min,max for all x values for needed variables
allx=np.unique(xCell[:])
N_mean = np.zeros(allx.shape)
N_min = np.zeros(allx.shape)
N_max = np.zeros(allx.shape)
q_mean = np.zeros(allx.shape)
q_min = np.zeros(allx.shape)
q_max = np.zeros(allx.shape)
for i in range(len(allx)):
    N_mean[i] = N[ xCell == allx[i] ].mean()
    N_min[i] = N[ xCell == allx[i] ].min()
    N_max[i] = N[ xCell == allx[i] ].max()
    q_mean[i] = q[ xCell == allx[i] ].mean()
    q_min[i] = q[ xCell == allx[i] ].min()
    q_max[i] = q[ xCell == allx[i] ].max()

print "start plotting."

fig = plt.figure(1, facecolor='w')
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312, sharex=ax1)
ax3 = fig.add_subplot(313, sharex=ax1)

if options.A5:
   G_x, G_Nmean, G_Nmin, G_Nmax, G_qmean, G_qmin, G_qmax, G_Qmax = np.loadtxt('tuning_A5.txt', skiprows=1, unpack=True, delimiter=",")

if options.A3:
   G_x, G_Nmean, G_Nmin, G_Nmax, G_qmean, G_qmin, G_qmax, G_Qmax = np.loadtxt('tuning_A5.txt', skiprows=1, unpack=True, delimiter=",")

if options.A3 or options.A5:
   # plot GLADS data
   lw = 3  # lineweight to use
   ax1.plot(G_x/1000.0, G_Nmin/1.0e6, 'g--', linewidth=lw)
   ax1.plot(G_x/1000.0, G_Nmean/1.0e6, 'g-', linewidth=lw, label='GLADS mean/range')
   ax1.plot(G_x/1000.0, G_Nmax/1.0e6, 'g--', linewidth=lw)

   ax2.plot(G_x/1000.0, G_qmin, 'g--', linewidth=lw)
   ax2.plot(G_x/1000.0, G_qmean, 'g-', linewidth=lw, label='GLADS mean/range')
   ax2.plot(G_x/1000.0, G_qmax, 'g--', linewidth=lw)

   ax3.plot(G_x/1000.0, G_Qmax, 'g--', linewidth=lw, label='GLADS max')

# panel 1: effective pressure
plt.sca(ax1)
#plt.plot(x, N[ind] / 1.0e6, '.-g')  # this just plots the centerline profile
plt.plot(allx/1000.0, N_mean / 1.0e6, '-b', label='MPAS mean/range')
plt.plot(allx/1000.0, N_min / 1.0e6, '--b')
plt.plot(allx/1000.0, N_max / 1.0e6, '--b')
plt.xlabel('X-position (km)')
plt.ylabel('effecive pressure (MPa)')
plt.xlim( (0, 100.0) )
plt.grid(True)
plt.legend(loc='best')

# panel 2: sheet flux
plt.sca(ax2)
#plt.plot(x, np.absolute(h[ind] * u[ind]), '.-g') # this plots centerline profile
plt.plot(allx/1000.0, np.absolute(q_mean), '-b', label='MPAS mean/range')
plt.plot(allx/1000.0, np.absolute(q_min), '--b')
plt.plot(allx/1000.0, np.absolute(q_max), '--b')
plt.xlabel('X-position (km)')
plt.ylabel('sheet water flux (m^2/s)')
plt.grid(True)
plt.legend(loc='best')

# panel 3: channel flux
plt.sca(ax3)
try:
   channelDischarge = f.variables['channelDischarge'][time_slice,:]
   allxEdge=np.unique(xEdge[:])
   Q_max = np.zeros(allxEdge.shape)
   Q_sum = np.zeros(allxEdge.shape)
   for i in range(len(allxEdge)):
     Q_max[i] = np.absolute(channelDischarge[ xEdge == allxEdge[i] ]).max()
     Q_sum[i] = np.absolute(channelDischarge[ xEdge == allxEdge[i] ]).sum()

   plt.plot(allxEdge/1000.0, np.absolute(Q_max), 'bx', label='MPAS max')
   plt.plot(allxEdge/1000.0, np.absolute(Q_sum), 'bo', label='MPAS sum')
except:
   print "Skipping plotting of channel output."

plt.xlabel('X-position (km)')
plt.ylabel('channel water flux (m^3/s)')
plt.grid(True)
plt.legend(loc='best')





# plot how close to SS we are
fig = plt.figure(2, facecolor='w')
ax1 = fig.add_subplot(211)
for i in ind:
    plt.plot(days/365.0, f.variables['waterThickness'][:,i])
plt.xlabel('Years since start')
plt.ylabel('water thickness (m)')
plt.grid(True)

ax = fig.add_subplot(212, sharex=ax1)
for i in ind:
    plt.plot(days/365.0, f.variables['effectivePressure'][:,i]/1.0e6)
plt.xlabel('Years since start')
plt.ylabel('effective pressure (MPa)')
plt.grid(True)



print "plotting complete"

plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

