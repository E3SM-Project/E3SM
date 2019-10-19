#!/usr/bin/env python

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

# test time evolution

plt.figure()

s = np.linspace(0,1,100)
x = s*25.0
y = -s*10.0

plt.plot(x, y, 'k-', lw=2, label='initial')


ds = xr.open_mfdataset('output.nc')

for atime in np.arange(len(ds.Time)):
    ds = ds.drop(np.setdiff1d(ds.variables.keys(), ['yCell','ssh','xtime']))
    ymean = ds.isel(Time=atime).groupby('yCell').mean(dim=xr.ALL_DIMS)
    x = ymean.yCell.values/1000.0
    y = ymean.ssh.values
    timelabel = str(ds.xtime.isel(Time=atime).values).strip()[11:]
    print('Doing {}'.format(timelabel))
    plt.plot(x, y, label=timelabel)

plt.ylim(-12,0)
plt.xlim(0,28)

plt.legend(fontsize=8)
plt.savefig('test_time_adv.png')

## test dependence on del4
#plt.figure()
#
#s = np.linspace(0,1,100)
#x = s*25.0
#y = -s*10.0
#
#plt.plot(x, y, 'k-', lw=2, label='initial')
#
#for strg in ['nodel4', '2e4', '1e6', '2e6']:
#    ds = xr.open_mfdataset('output' + strg + '.nc')
#    ds = ds.drop(np.setdiff1d(ds.variables.keys(), ['yCell','ssh']))
#    ymean = ds.isel(Time=-1).groupby('yCell').mean(dim=xr.ALL_DIMS)
#    x = ymean.yCell.values/1000.0
#    y = ymean.ssh.values
#    plt.plot(x, y, label=strg)
#
#plt.legend()
#plt.savefig('test_del4.png')
#
#
