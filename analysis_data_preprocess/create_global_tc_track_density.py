#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The observational TC data are obtained from Prof. Kerry Emanuel’s website (https://emanuel.mit.edu/products). Click “Global Tropical Cyclone Data in NETCDF format (updated through 2018)” the files will be downloaded. The TC locations are tracked 6 hourly.

Another data source is IBTrACS from NOAA, which has 3hourly track data.

Note, The tc density maps created from both data source have identical distribution, but values are twice from IBTrACS compared to MIT due to higher time frequency. For now the TC analysis has been using 6hourly E3SM output, therefore the MIT data is being used for consistency.
"""

import os
import sys
import warnings
from datetime import datetime, timedelta

if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
from netCDF4 import Dataset as netcdffile

all_lon = []
all_lat = []

###################################################
origin_path = '/Users/zhang40/Documents/ACME/e3sm_tc_diags'
start_yr = 1979
end_yr = 2018

ibtracs = False #MIT data
ibtracs = True


if ibtracs:
    data_name = 'IBTrACS'
    print('Use IBTrACS 3hourly')
    basins = ["NA", "WP", "EP", "NI", "SI", "SP"]
    for basin in basins:
        nc = netcdffile(
            os.path.join(origin_path, "IBTrACS.{}.v04r00.nc".format(basin))
        )
        latmc = np.squeeze(nc["lat"][:, :])
        longmc = np.squeeze(nc["lon"][:, :])
        time = np.squeeze(nc["time"][:, :])
        yearic = np.zeros((np.shape(latmc)[0]))

        for i in range(0, np.shape(time)[0]):
            for j in range(0, np.shape(time)[1]):
                if time[i, j] is not np.ma.masked:
                    day_hurr = datetime(1858, 11, 17, 0, 0, 0) + timedelta(time[i, j])
                    yearic[i] = day_hurr.year
                    if yearic[i] >= start_yr and yearic[i] <= end_yr and np.abs(latmc[i,j]) > 0 and np.abs(longmc[i,j]) > 0:
                        if(longmc[i, j]<0):
                            longmc[i, j] = 360 + longmc[i, j]
                        all_lat.append(latmc[i, j])
                        all_lon.append(longmc[i, j])
    print(len(all_lat))


else:
    print('MIT 6hrly')
    data_name = 'MIT'
    basins = ['at', 'ep','wp','io','sh']
    for basin in basins:
        nc = netcdffile('{}/{}tracks.nc'.format(origin_path, basin))
        latmc = np.squeeze(nc['latmc'][:,:])
        longmc = np.squeeze(nc['longmc'][:,:])
        vsmc = np.squeeze(nc['vsmc'][:,:])
        yearic = np.squeeze(nc['yearic'][:])

        for i in range(0,np.shape(latmc)[0]):
            for j in range(0,np.shape(latmc)[1]):

                if yearic[j] >=start_yr and yearic[j]<=end_yr and np.abs(latmc[i,j]) > 0 and np.abs(longmc[i,j]) > 0:
                   all_lon.append(longmc[i,j])
                   all_lat.append(latmc[i,j])
    print(len(all_lon))

all_lat = np.asarray(all_lat)
all_lon = np.asarray(all_lon)

all_lat_new = np.round(all_lat,2)
all_lon_new = np.round(all_lon,2)

import pandas as pd

raw_data = {
        'lon': all_lon_new,
        'lat': all_lat_new}

df = pd.DataFrame.from_dict(raw_data)
print(data_name)
out_file = '/Users/zhang40/Documents/ACME/e3sm_tc_diags/cyclones_hist_{}_{}_{}.csv'.format(data_name,start_yr, end_yr)
with open(out_file, 'w') as file:
    file.write('start	{}\n'.format(len(all_lon)))
    df.to_csv(file, header=False, sep='\t')

###################################################
# Convert the .csv file to .nc by calling tempest-extremes from command line:
#tempestextremes/bin/HistogramNodes --in cyclones_hist_IBTrACS_1979_2018.csv --iloncol 2 --ilatcol 3 --out cyclones_hist_IBTrACS_1979_2018.nc
