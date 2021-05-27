#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 11:06:41 2019

@author: bala635
The observational TC data are obtained from Prof. Kerry Emanuel’s website (https://emanuel.mit.edu/products). Click “Global Tropical Cyclone Data in NETCDF format (updated through 2018)” the files will be downloaded.
"""

import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
from netCDF4 import Dataset as netcdffile

###################################################
###################################################

all_lon = []
all_lat = []

###################################################

nc = netcdffile('./attracks.nc')
latmc = np.squeeze(nc['latmc'][:,:])
longmc = np.squeeze(nc['longmc'][:,:])
vsmc = np.squeeze(nc['vsmc'][:,:])
yearic = np.squeeze(nc['yearic'][:])

for i in range(0,np.shape(latmc)[0]):
    for j in range(0,np.shape(latmc)[1]):

        if yearic[j] > 1980 and np.abs(latmc[i,j]) > 0 and np.abs(longmc[i,j]) > 25:

           all_lon.append(longmc[i,j])
           all_lat.append(latmc[i,j])

###################################################

nc = netcdffile('./eptracks.nc')
latmc = np.squeeze(nc['latmc'][:,:])
longmc = np.squeeze(nc['longmc'][:,:])
vsmc = np.squeeze(nc['vsmc'][:,:])
yearic = np.squeeze(nc['yearic'][:])

for i in range(0,np.shape(latmc)[0]):
    for j in range(0,np.shape(latmc)[1]):

        if yearic[j] > 1980 and np.abs(latmc[i,j]) > 0 and np.abs(longmc[i,j]) > 25:

           all_lon.append(longmc[i,j])
           all_lat.append(latmc[i,j])

###################################################

nc = netcdffile('./wptracks.nc')
latmc = np.squeeze(nc['latmc'][:,:])
longmc = np.squeeze(nc['longmc'][:,:])
vsmc = np.squeeze(nc['vsmc'][:,:])
yearic = np.squeeze(nc['yearic'][:])

for i in range(0,np.shape(latmc)[0]):
    for j in range(0,np.shape(latmc)[1]):

        if yearic[j] > 1980 and np.abs(latmc[i,j]) > 0 and np.abs(longmc[i,j]) > 25:

           all_lon.append(longmc[i,j])
           all_lat.append(latmc[i,j])

###################################################

nc = netcdffile('./iotracks.nc')
latmc = np.squeeze(nc['latmc'][:,:])
longmc = np.squeeze(nc['longmc'][:,:])
vsmc = np.squeeze(nc['vsmc'][:,:])
yearic = np.squeeze(nc['yearic'][:])

for i in range(0,np.shape(latmc)[0]):
    for j in range(0,np.shape(latmc)[1]):

        if yearic[j] > 1980 and np.abs(latmc[i,j]) > 0 and np.abs(longmc[i,j]) > 25:

           all_lon.append(longmc[i,j])
           all_lat.append(latmc[i,j])

###################################################

nc = netcdffile('./shtracks.nc')
latmc = np.squeeze(nc['latmc'][:,:])
longmc = np.squeeze(nc['longmc'][:,:])
vsmc = np.squeeze(nc['vsmc'][:,:])
yearic = np.squeeze(nc['yearic'][:])

for i in range(0,np.shape(latmc)[0]):
    for j in range(0,np.shape(latmc)[1]):

        if yearic[j] > 1980 and np.abs(latmc[i,j]) > 0 and np.abs(longmc[i,j]) > 25:

           all_lon.append(longmc[i,j])
           all_lat.append(latmc[i,j])

###################################################

all_lat = np.asarray(all_lat)
all_lon = np.asarray(all_lon)

all_lat_new = np.round(all_lat,2)
all_lon_new = np.round(all_lon,2)

import pandas as pd

raw_data = {'lon': all_lon_new,
        'lat': all_lat_new}

df = pd.DataFrame.from_dict(raw_data)
# columns = ['SST', 'STRAT','TCHP','TDY_SALT','TDY_NOSALT','INT','LAT','LON','SHR','DIV','PER','DELTA_INT']
df.to_csv('cyclones_all_obs.csv', sep='\t')

###################################################
# Convert the .csv file to .nc by calling tempest-extremes from command line:
#tempestextremes/bin/HistogramNodes --in cyclones_all_obs.csv --iloncol 3 --ilatcol 4 --out cyclones_all_obs.nc
