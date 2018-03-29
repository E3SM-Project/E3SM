#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np

data = Dataset('./ais20km_bmb_background_anomaly.nc','r+')

for i in range(100):

    time_string = data.variables['xtime'][i,:]

    time_string[0:4] = np.array(list(str(i+1).zfill(4)))

    print time_string

    data.variables['xtime'][i,:] = time_string


data.close()
