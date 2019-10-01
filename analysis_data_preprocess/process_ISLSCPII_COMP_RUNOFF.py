#This is a data conversion script to convert ISLSCPII_COMP_RUNOFF from ASCII to netcdf files
#for use with e3sm_diags package
#The data was taken from the ASCII format original data available from NASA EARTH DATA
#https://daac.ornl.gov/ISLSCP_II/guides/comp_runoff_monthly_xdeg.html
# Aug 1st, 2019 by Jill Chengzhu Zhang

import cdms2
import MV2
import glob
import numpy
import cdutil
from calendar import monthrange


basedir = '/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ISLSCPII_GRDC/'
datapath = basedir  + 'original_data/data/comp_runoff_hdeg/'
file_name = 'comp_runoff_hd'

start_time = '198601'
end_time = '199512'
#Missing data cells are given the value of -99 over water, and -88 over land 
#(mainly Greenland and Antarctica).
missing = -99
missing2 = -88

nlat = 360
nlon = 720
#nlat = 180
#nlon =360

#Create axes, and a variable to hold 12 months of data
lat = cdms2.createAxis([-89.25 + 0.5*i for i in range(nlat)],id='lat')
lon = cdms2.createAxis([-179.75 + 0.5*i for i in range(nlon)],id='lon')
out_file = cdms2.open(basedir + 'time_series/QRUNOFF_' + start_time + '_' + end_time + '.nc','w')
start = "1986-01-01"
nTimes = 12 * 10
time = cdms2.createAxis(range(nTimes))
time.units = "months since {}".format(start)
time.designateTime()
time.id = "time"
data = MV2.zeros([nTimes, nlat, nlon])

print(data.shape)

count = 0
for year in range(1986,1996):
    for month in range(12):
        print(year,month)
        days_per_month = monthrange(year, month+1)[1]
        print(days_per_month)
        imon = '{0:02d}'.format(month+1)
        fname = datapath + file_name +'_' + str(year)+imon+'.asc'
        data[count,:,:] = numpy.loadtxt(fname, skiprows=6) / days_per_month * 24. * 3600.  #converted from mm/month to kg/m2/s
        count = count+1

print(count)
#data = MV2.masked_where(data == missing, data)
#data = MV2.masked_where(data == missing2, data)
data = MV2.masked_where(data < 0 , data)
out_data = cdms2.createVariable(data, axes =[time,lat,lon], id = "mrro")
#out_data.units = "mm/month"   
out_data.units = "kg/m2/s"   #converted from mm/month to kg/m2/s
out_data.id = "QRUNOFF" #runoff flux for cmip
cdutil.setTimeBoundsMonthly(out_data)
out_file.dataname = 'ISLSCP-GRDC'
out_file.long_dataname = "ISLSCP II UNH/GRDC Composite Monthly Runoff"
out_file.comment = "The files in this data set comprise river discharge point data from the Global Runoff Data Centre (GRDC) for all 390 stations with measurements in the period 1986-1995 and 2) global, gridded, monthly composite runoff time series for the period 1986-1995, produced at UNH using the GRDC river discharge data. The original data was provided in ASCII format, it is converted to netcdf for using in e3sm_diags"
out_file.years = "1986-1995"
out_file.write(out_data)



#import cdutil
#seasonal = cdutil.DJF.climatology(out_file)
#print(seasonal)


