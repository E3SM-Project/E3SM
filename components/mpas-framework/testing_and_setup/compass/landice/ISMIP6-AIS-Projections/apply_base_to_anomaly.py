#!/usr/bin/env python

# script to apply anomaly SMB field to our base SMB

from __future__ import absolute_import, division, print_function, unicode_literals
import netCDF4
import numpy as np

f = netCDF4.Dataset('SMB_NorESM_2.6_1995-2100.nc','r+')
nt = len(f.dimensions['Time'])
smb = f.variables['sfcMassBal']

fbase = netCDF4.Dataset('../test/ais2km_100yr_spinup.nc','r')
smbBase = fbase.variables['sfcMassBal'][0,:]

for i in range(nt):
    print(i)

    smb[i,:] = smb[i,:] + smbBase

f.close()
fbase.close()

