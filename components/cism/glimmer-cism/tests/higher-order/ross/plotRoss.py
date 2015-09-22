#!/usr/bin/env python
# This script plots the results of the Ross Ice Shelf Experiments.
# The plots compare the model output with measured velocities from the 
# Ross Ice Shelf Geophysical and Glaciological Survey 1973-78 (RIGGS).
# Before running this script run runGlimmer.py to generate the model output.
# See the accompanying README file for more information.
# Written April 5, 2010 by Glen Granzow at the University of Montana.

import os
import numpy
from optparse import OptionParser
from matplotlib import pyplot, colors
from netCDF import *

########## PARSE ANY COMMAND LINE ARGUMENTS ##########

parser = OptionParser()
parser.add_option('-c','--cmap',dest='colormap',type='string',default='jet',help='Specify a color map')
parser.add_option('-m','--mask',dest='use_mask',action='store_true',default=False,help='Mask the grounded region')
parser.add_option('-n','--ncontours',dest='ncontours',type='int',default=0,help='Plot n filled contours instead of a continuous color plot')
parser.add_option('-v','--vmax',dest='vmax',type='float',default=0,help='Set upper limit of color scale')
options, args = parser.parse_args()
print 'options:', options

# cmap_discretize returns a discretized version of a colormap
# From http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations 
def cmap_discretize(cmap, N):
     from scipy import interpolate
     cdict = cmap._segmentdata.copy()
     colors_i = numpy.linspace(0,1,N)
     indices  = numpy.linspace(0,1,N+1)
     for key in ('red','green','blue'):
         D = numpy.array(cdict[key])
         I = interpolate.interp1d(D[:,0], D[:,1])
         clrs = I(colors_i)
         A = numpy.zeros((N+1,3), float)
         A[:,0] = indices
         A[1:,1] = clrs
         A[:-1,2] = clrs
         L = []
         for l in A:
             L.append(tuple(l))
         cdict[key] = tuple(L)
     return colors.LinearSegmentedColormap('colormap',cdict,1024)

########## PART I: READ THE INPUT FILES ##########

# Read the Glimmer output file
filename = os.path.join('output','ross.out.nc')
inputfile1 = NetCDFFile(filename,'r')
velnorm = numpy.array(inputfile1.variables['velnorm'][0,0,2:-2,2:-2])
inputfile1.close()

if options.use_mask:
  filename = os.path.join('output','ross.nc')
  inputfile2 = NetCDFFile(filename,'r')
  mask = numpy.array(inputfile2.variables['kinbcmask'][0,2:-2,2:-2])
  inputfile2.close()

# Read the grid coordinates 
filename = os.path.join('data','111by147Grid.dat')
inputfile3 = open(filename)
nx,ny = 147,111
x,y = list(),list()
for i in range(4): 
  inputfile3.readline()
for i in range(ny):
  y.append(float(inputfile3.readline()))
for i in range(3):
  inputfile3.readline()
for i in range(nx):
  x.append(float(inputfile3.readline()))
inputfile3.close()

# Read the RIGGS data
filename = os.path.join('data','riggs_clean.dat')
inputfile4 = open(filename)
latitude,longitude,riggs_velocity,glimmer_velocity = list(),list(),list(),list()
lat0,lon0,vel0 = list(),list(),list()
for line in inputfile4:
  tokens = [float(data) for data in line.split()]
# Convert from degrees, minutes and seconds
  lat = -(tokens[3] + tokens[4]/60 + tokens[5]/60**2)
  lon = -(tokens[6] + tokens[7]/60 + tokens[8]/60**2)*tokens[9]
  if y[0] < lat < y[-1] and x[0] < lon < x[-1]:
#   Interpolate the Glimmer data onto this point
    for i in range(nx):
      if lon < x[i+1]: break
    for j in range(ny):
      if lat < y[j+1]: break
    alpha = (lon-x[i])/(x[i+1]-x[i])
    beta  = (lat-y[j])/(y[j+1]-y[j])
    v = (1-alpha)*((1-beta)*velnorm[j,i]  +beta*velnorm[j+1,i]) \
       +   alpha *((1-beta)*velnorm[j,i+1]+beta*velnorm[j+1,i+1])
    if v != 0:
      latitude.append(lat)
      longitude.append(lon)
      riggs_velocity.append(tokens[10])
      glimmer_velocity.append(v)
    else:
      print 'Glimmer velocity is zero at (%g, %g); RIGGS velocity is %g' % (lon,lat,tokens[10])
      lat0.append(lat)
      lon0.append(lon)
      vel0.append(tokens[10])
  else:
    print 'RIGGS data is off the grid at (%g, %g); RIGGS velocity is %g' % (lon,lat,tokens[10])
print 'These points are not included in Chi-squared.'
inputfile4.close()

########## PART II: PRESENT THE OUTPUT ##########

sigma = 30
X2 = sum([((v1-v2)/sigma)**2 for v1,v2 in zip(riggs_velocity,glimmer_velocity)])
print
print 'Chi-squared for',len(riggs_velocity),'points is', X2
print 'The maximum velocity from Glimmer/CISM is', numpy.max(velnorm)
print 'The maximum velocity from the RIGGS data is', max(riggs_velocity)

# Create a scatter plot
pyplot.figure(1)
pyplot.clf()
pyplot.plot([0,1200],[0,1200],color='red',linewidth=2,alpha=0.5)
pyplot.scatter(riggs_velocity,glimmer_velocity,color='blue')
pyplot.axis([0,1050,0,1400])
pyplot.xlabel('Measured velocity from RIGGS (meters/year)')
pyplot.ylabel('Model velocity from Glimmer/CISM (meters/year)')
pyplot.title('Ross Ice Shelf Experiment')

# Create a color plot of Glimmer and RIGGS velocities
pyplot.figure(2)
pyplot.clf()
plot_discarded_points = False
# The dimensions of X and Y for pcolor should be one greater than those of the data
midx = sum(map(numpy.array,(x[:-1],x[1:])))/2
midy = sum(map(numpy.array,(y[:-1],y[1:])))/2
x = x[:1] + list(midx) + x[-1:]
y = y[:1] + list(midy) + y[-1:]
x,y = numpy.meshgrid(x,y)
# Set plot parameters
vmax = numpy.max(velnorm)
if options.vmax > 0:
  vmax = options.vmax
norm  = colors.Normalize(0,vmax)
cmap  = eval('pyplot.cm.' + options.colormap)
ticks = None
if options.ncontours > 1:
  cmap = cmap_discretize(cmap, options.ncontours)
  ticks = numpy.arange(options.ncontours+1) * vmax/options.ncontours
if options.use_mask:
  velnorm = numpy.ma.masked_array(velnorm,mask==1)
# Make the actual plot
pyplot.pcolor(x,y,velnorm,norm=norm,cmap=cmap)
if plot_discarded_points:
  pyplot.scatter(lon0,lat0,s=50,c='white')
pyplot.scatter(longitude,latitude,s=50,c=riggs_velocity,norm=norm,cmap=cmap)
pyplot.axis([x[0,0],x[-1,-1],y[0,0],y[-1,-1]])
pyplot.axes().set_aspect('equal')
pyplot.colorbar(orientation='vertical',fraction=0.05,pad=0.01,shrink=0.93,ticks=ticks)
pyplot.title('Ross Ice Shelf Experiment - Velocity (meters/year)')
print
print 'FIGURE 2:'
print 'If Glimmer/CISM perfectly predicted the velocities measured by the Ross'
print 'Ice Shelf Geophysical and Glaciological Survey (RIGGS), the colors in'
print 'the circles (which represent the RIGGS velocities) would be the same as'
print 'the color of the background (which represents the velocities calculated'
print 'by Glimmer/CISM).'
print
pyplot.show()
