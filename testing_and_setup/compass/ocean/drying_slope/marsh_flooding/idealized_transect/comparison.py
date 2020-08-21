#!/usr/bin/env python
# coding: utf-8
# Authors: Zhendong Cao, Phillip Wolfram
# Date: May 2020

# import libs
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
#import cv2
# plot tidal amplitude at open boundary

fname1 = 'NoVegetation'
fname2 = 'ConstantVegManning'
fname3 = 'VegManningEquation'

# mpas-o
ds1 = xr.open_dataset(fname1+'.nc')
ds2 = xr.open_dataset(fname2+'.nc')
ds3 = xr.open_dataset(fname3+'.nc')

x = np.linspace(0, 1.0, len(ds1.xtime))*48.0
ympas1 = ds1.ssh.where(ds1.tidalInputMask).mean('nCells').values
ympas2 = ds2.ssh.where(ds2.tidalInputMask).mean('nCells').values
ympas3 = ds3.ssh.where(ds3.tidalInputMask).mean('nCells').values

plt.plot(x, ympas1, marker='o', label=fname1)
plt.plot(x, ympas2, marker='^', label=fname2)
plt.plot(x, ympas3, marker='*', label=fname3)

# analytical
x = np.linspace(0,48.0,100)
y = 1.5*np.sin(x*np.pi/12.0) - 3.0
plt.plot(x, y, 'k-', lw=3, label='analytical')

plt.legend(frameon=False, loc='upper right', fontsize=8)
plt.ylabel('Tidal amplitude (m)')
plt.xlabel('Time (hrs)')
plt.suptitle('Tidal amplitude forcing for test cases and analytical')
plt.savefig('tidalcomparison.png',dpi=300)
plt.close()

# extract bottomDepth
bottomDepth = ds1.bottomDepth.values.reshape((116, 4))[:,1]
x_along = np.linspace(0,8,116)
# plot water elevation at selected time slices
times = np.linspace(0.0,2.0,41)
ii = 0
figname_prefix = 'Waterlevel'
for atime in times:
    plt.figure(figsize=(12,6))
    plt.xlim(0, 8)
    plt.ylim(-4, 0)
    plottime = int((float(atime)/0.2)*24.0)
    timestr = ds1.isel(Time=plottime).xtime.values
    ymean = ds1.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
    x = ymean.yCell.values/1000.0
    y = ymean.ssh.values
    plt.plot(x, y,'-k',label=fname1)

    ymean = ds2.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
    x = ymean.yCell.values/1000.0
    y = ymean.ssh.values
    plt.plot(x, y,'-b',label=fname2)

    ymean = ds3.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
    x = ymean.yCell.values/1000.0
    y = ymean.ssh.values
    plt.plot(x, y,'-g',label=fname3)

    plt.fill_between(x=x_along,y1=-10*np.ones(x_along.shape), y2=-bottomDepth, color='grey',zorder=9)
    plt.xlabel('Cross shore distance (km)')
    plt.ylabel('Bathymetry (m)')
    atime='%.2f'%atime
    plt.title('time='+atime+'days')
    plt.legend(frameon=False, loc='upper right')
    figname = '{}{:02d}.png'.format(figname_prefix,ii)
    plt.savefig(figname, dpi=100)
    plt.close()
    print('Time '+atime+ ' done!')
    ii += 1
# make a movie
os.system("ffmpeg -r 1 -i Waterlevel%02d.png -vcodec mpeg4 -y movie.mp4")
