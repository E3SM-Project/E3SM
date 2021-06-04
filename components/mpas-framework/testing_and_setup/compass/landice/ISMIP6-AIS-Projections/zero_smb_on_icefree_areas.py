#!/usr/bin/env python

# script to apply anomaly SMB field to our base SMB

from __future__ import absolute_import, division, print_function, unicode_literals
import netCDF4
import numpy as np

fic = netCDF4.Dataset('/global/cscratch1/sd/hoffman2/ISMIP6/expt5-pio1/ais2km_100yr_spinup.nc','r')
thk0 = fic.variables['thickness'][0,:]
bed = fic.variables['bedTopography'][0,:]
mask = (thk0>0)  # only keep SMB where ice was originally
fic.close()

f = netCDF4.Dataset('SMB_NorESM_2.6_1995-2100.nc','r+')
nt = len(f.dimensions['Time'])
smb = f.variables['sfcMassBal']

for i in range(nt):
    print(i)

    smb[i,:] = smb[i,:] * mask

f.close()

