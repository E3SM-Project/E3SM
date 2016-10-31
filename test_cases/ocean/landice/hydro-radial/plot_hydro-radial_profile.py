#!/usr/bin/env python
'''
Plots profiles for hydro-margin test case
'''
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

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

if not options.time:
	print "No time provided. Using time -1."
        time_slice = -1
else:
        time_slice = int(options.time)



f = netCDF4.Dataset(options.filename,'r')
#xtime = f.variables['xtime'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
xEdge = f.variables['xEdge'][:]
yEdge = f.variables['yEdge'][:]
h = f.variables['waterThickness'][time_slice,:]
u = f.variables['waterVelocityCellX'][time_slice,:]
P = f.variables['waterPressure'][time_slice,:]
N = f.variables['effectivePressure'][time_slice,:]
div = f.variables['divergence'][time_slice,:]
#H = f.variables['thickness'][time_slice,:]
opening = f.variables['openingRate'][time_slice,:]
closing = f.variables['closingRate'][time_slice,:]
melt = f.variables['basalMeltInput'][time_slice,:]
sliding = f.variables['basalSpeed'][time_slice,:]
days = f.variables['daysSinceStart'][:]
xtime = f.variables['xtime'][:]

print "Total number of time levels=", len(days)
print "Using time slice", time_slice, " which is year ", days[time_slice]/365.0
print "xtime=", ''.join(xtime[time_slice,:])

print "Attempting to read thickness field from landice_grid.nc."
fin = netCDF4.Dataset("landice_grid.nc",'r')
H = fin.variables['thickness'][0,:]


# Find center row  - currently files are set up to have central row at y=0
unique_ys=np.unique(yCell[:])
centerY=unique_ys[len(unique_ys)/2]
print "number of ys, center y index, center Y value", len(unique_ys), len(unique_ys)/2, centerY
ind = np.nonzero(yCell[:] == centerY)
x = xCell[ind]/1000.0

print "start plotting."

fig = plt.figure(1, facecolor='w')
# water thickness
ax1 = fig.add_subplot(121)
#plt.plot(x, H[ind]*917.0*9.81/1.0e5, '.-')
plt.plot(x, h[ind], '.-')
plt.xlabel('X-position (km)')
plt.ylabel('water depth (m)')
plt.grid(True)

# water pressure
ax = fig.add_subplot(122, sharex=ax1)
plt.plot(x, H[ind]*910.0*9.80616 / 1.0e5, '.-')
plt.plot(x, P[ind] / 1.0e5, '.--')
plt.xlabel('X-position (km)')
plt.ylabel('water pressure (bar)')
plt.grid(True)


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

# plot opening/closing rates
fig = plt.figure(3, facecolor='w')

nplt=5

ax = fig.add_subplot(nplt,1,1)
plt.plot(x, opening[ind], 'r', label='opening')
plt.plot(x, closing[ind], 'b', label='closing')
plt.plot(x, melt[ind] / 1000.0, 'g', label='melt')
plt.xlabel('X-position (km)')
plt.ylabel('rate (m/s)')
plt.legend()
plt.grid(True)

# SS N=f(h0
ax = fig.add_subplot(nplt,1,2)
plt.plot(x, N[ind]/1.0e6, '.-', label='modeled transient to SS')
plt.plot(x, (opening[ind]/(0.04*3.1709792e-24*h[ind]))**0.3333333 / 1.0e6, '.--r', label='SS N=f(h)')  # steady state N=f(h) from the cavity evolution eqn
plt.xlabel('X-position (km)')
plt.ylabel('effective pressure (MPa)')
#plt.ylim((0.0, 0.1))
plt.grid(True)
plt.legend()

ax = fig.add_subplot(nplt,1,3)
plt.plot(x, u[ind])
plt.ylabel('water velocity (m/s)')
plt.grid(True)


ax = fig.add_subplot(nplt,1,4)
plt.plot(x, u[ind]*h[ind])
plt.ylabel('water flux (m2/s)')
plt.grid(True)

ax = fig.add_subplot(nplt,1,5)
plt.plot(x, div[ind])
plt.plot(x, melt[ind] / 1000.0, 'g', label='melt')
plt.ylabel('divergence (m/s)')
plt.grid(True)

# optional - check velo field correctness
#fig = plt.figure(4, facecolor='w')
#plt.plot(x, sliding[ind])
#plt.grid(True)


# plot some edge quantities
inde = np.nonzero(yEdge[:] == centerY)
xe = xEdge[inde]/1000.0
ve = f.variables['waterVelocity'][time_slice,:]
#k = f.variables['effectiveConducEdge'][time_slice,:]
dphie = f.variables['hydropotentialBaseSlopeNormal'][time_slice,:]
he = f.variables['waterThicknessEdgeUpwind'][time_slice,:]
fluxe = f.variables['waterFluxAdvec'][time_slice,:]

fig = plt.figure(5, facecolor='w')
nplt=5

ax1 = fig.add_subplot(nplt,1,1)
plt.plot(xe, dphie[inde],'.')
plt.ylabel('dphidx edge)')
plt.grid(True)

ax = fig.add_subplot(nplt,1,2, sharex=ax1)
plt.plot(x, P[ind],'x')
plt.ylabel('dphidx edge)')
plt.grid(True)

ax = fig.add_subplot(nplt,1,3, sharex=ax1)
plt.plot(xe, ve[inde],'.')
plt.ylabel('vel edge)')
plt.grid(True)

ax = fig.add_subplot(nplt,1,4, sharex=ax1)
plt.plot(xe, he[inde],'.')
plt.plot(x, h[ind],'x')
plt.ylabel('h edge)')
plt.grid(True)

ax = fig.add_subplot(nplt,1,5, sharex=ax1)
plt.plot(xe, fluxe[inde],'.')
plt.ylabel('flux edge)')
plt.grid(True)


# ==========
# Make plot similar to Bueler and van Pelt Fig. 5

# get thickness/pressure at time 0 - this should be the nearly-exact solution interpolated onto the MPAS mesh
h0 = f.variables['waterThickness'][0,:]
P0 = f.variables['waterPressure'][0,:]
hasice = (sliding>0.0)  # assuming sliding has been zeroed where there is no ice, so we don't need to get the thickness field

Werr = np.absolute(h - h0)
Perr = np.absolute(P - P0)
dcEdge= f.variables['dcEdge'][:]
dx = dcEdge.mean()  # ideally should restrict this to edges with ice

fig = plt.figure(6, facecolor='w')

ax = fig.add_subplot(2,1,1)
plt.plot(dx, Werr[hasice].mean(), 's', label='avg W err')
plt.plot(dx, Werr[hasice].max(), 'x', label='max W err')
ax.set_yscale('log')
plt.grid(True)
plt.legend()
plt.xlabel('delta x (m)')
plt.ylabel('error in W (m)')

ax = fig.add_subplot(2,1,2)
plt.plot(dx, Perr[hasice].mean()/1.0e5, 's', label='avg P err')
plt.plot(dx, Perr[hasice].max()/1.0e5, 'x', label='max P err')
ax.set_yscale('log')
plt.grid(True)
plt.legend()
plt.xlabel('delta x (m)')
plt.ylabel('error in P (bar)')


print "plotting complete"
plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

