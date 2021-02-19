import cdms2
from pathlib import Path
#from acme_diags.driver.utils import climo
from acme_diags.driver.utils import general
import climo
import diurnal_cycle
import numpy
import os

value = 0
cdms2.setNetcdfShuffleFlag(value) ## where value is either 0 or 1
cdms2.setNetcdfDeflateFlag(value) ## where value is either 0 or 1
cdms2.setNetcdfDeflateLevelFlag(value) ## where value is a integer between 0 and 9 included

# In this example 3 hourly output at ARM sites are saved on h1 tape using namelist as follows:
#fincl2 = 'PS', 'Q', 'T', 'Z3', 'CLOUD', 'CONCLD', 'CLDICE', 'CLDLIQ', 'LS_FLXPRC', 'LS_FLXSNW', 'ZMFLXPRC', 'ZMFLXSNW', 'FREQR', 'REI', 'REL', 'CV_REFFICE', 'CV_REFFLIQ', 'LS_REFFRAIN', 'LS_REFFSNOW', 'PRECT', 'TMQ', 'PRECC', 'TREFHT', 'QREFHT', 'OMEGA','CLDTOT', 'LHFLX', 'SHFLX', 'FLDS', 'FSDS', 'FLNS', 'FSNS', 'FLNSC', 'FSDSC', 'FSNSC', 'AODVIS', 'AODABS'
#fincl2lonlat = '262.5e_36.6n','204.6e_71.3n','147.4e_2.0s','166.9e_0.5s','130.9e_12.4s','331.97e_39.09n'
data_path = '/Users/zhang40/Documents/ACME/e3sm_arm_diags/data/20200922.F2010SC5.ne30pg2_r05.armsites/'
out_path = '/Users/zhang40/Documents/ACME/e3sm_arm_diags/data/post-processed/20200922.F2010SC5.ne30pg2_r05.armsites/'
p = Path(data_path)
cmd = 'ncrcat -h '+data_path+'*h1*nc ' +data_path+'armsites_all_time.nc' 
os.popen(cmd).readlines()

filename = data_path+'armsites_all_time.nc'
fin = cdms2.open(filename)
variables = ['CLOUD', 'CONCLD', 'CLDICE', 'CLDLIQ', 'PRECT', 'TMQ', 'TREFHT', 'QREFHT', 'CLDTOT', 'LHFLX', 'SHFLX', 'FLDS', 'FSDS', 'FLNS', 'FSNS', 'FLNSC', 'FSDSC', 'FSNSC', 'AODVIS', 'AODABS', 'PS']
sites_info = {'sgp':[262.5,36.6], 'nsa':[204.6, 71.3], 'twpc1':[147.4,-2.], 'twpc2':[166.9, -0.5], 'twpc3':[130.9,-12.4], 'ena':[331.97, 39.09]}
sites = ['sgp','nsa','twpc1','twpc2','twpc3']

time_range = '000101_000112'


for site in sites:
    if sites_info[site][1] >0: 
        lon_lat = str(sites_info[site][0])+'e_'+str(sites_info[site][1])+'n'
    else:
        lon_lat = str(sites_info[site][0])+'e_'+str(abs(sites_info[site][1]))+'s'
    
    for variable in variables:
        fout_1 = cdms2.open(out_path+ variable+'_'+site+'_'+time_range+'.nc','w')
        var_name = variable + '_' +lon_lat
        
        print(var_name)
        var = fin(var_name,squeeze=1)
        var_time = var.getTime()
    
        if ' 0000-' in var_time.units:
            units = var_time.units
            fakeUnits = units.replace(' 0000-', ' 0001-')   # valid units for cdtime https://github.com/CDAT/cdms/issues/334
            var_time.units = fakeUnits
        var.id = variable
        lat = fin('lat_'+lon_lat)
        lat = lat[:]
        lat = cdms2.createVariable(lat)
        lat.id = 'lat'
        lat.long_name = "latitude"
        lat.units = "degrees_north"
        lon = fin('lon_'+lon_lat)
        lon = lon[:]
        lon = cdms2.createVariable(lon)
        lon.id = 'lon'
        lon.long_name = "longitude"
        lon.units = "degrees_east"
        fout_1.write(var)
        fout_1.write(lat)
        fout_1.write(lon)
        if var.getLevel():
            var1 = fin('PS'+ '_' +lon_lat,squeeze=1)  
            var1.id = 'PS'
            var2 = fin('P0')  
            var2.id = 'P0'
            var3 = fin('hyam')  
            var3.id = 'hyam'
            var4 = fin('hybm')  
            var4.id = 'hybm'
            fout_1.write(var1)
            fout_1.write(var2)
            fout_1.write(var3)
            fout_1.write(var4)
        fout_1.close()

        
