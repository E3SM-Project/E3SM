import os
from pathlib import Path

import cdms2

# subprocess.run('source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh', shell=True)

value = 0
cdms2.setNetcdfShuffleFlag(value)  # where value is either 0 or 1
cdms2.setNetcdfDeflateFlag(value)  # where value is either 0 or 1
cdms2.setNetcdfDeflateLevelFlag(
    value
)  # where value is a integer between 0 and 9 included

# A script to convert high frequency single point E3SM output to per-variable per-site netcdf files as input for ARM diagostics.
# In this example 3 hourly output at ARM sites are saved on h4 tape using namelist as follows:
# fincl2 = 'PS', 'Q', 'T', 'Z3', 'CLOUD', 'CONCLD', 'CLDICE', 'CLDLIQ', 'LS_FLXPRC', 'LS_FLXSNW', 'ZMFLXPRC', 'ZMFLXSNW', 'FREQR', 'REI', 'REL', 'CV_REFFICE', 'CV_REFFLIQ', 'LS_REFFRAIN', 'LS_REFFSNOW', 'PRECT', 'TMQ', 'PRECC', 'TREFHT', 'QREFHT', 'OMEGA','CLDTOT', 'LHFLX', 'SHFLX', 'FLDS', 'FSDS', 'FLNS', 'FSNS', 'FLNSC', 'FSDSC', 'FSNSC', 'AODVIS', 'AODABS'
# fincl2lonlat = '262.5e_36.6n','204.6e_71.3n','147.4e_2.0s','166.9e_0.5s','130.9e_12.4s','331.97e_39.09n'

# data_path = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210719.PhaseII.F20TR-P3.NGD.ne30pg2.compy/h4/"
# out_path = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/20210719.PhaseII.F20TR-P3.NGD.ne30pg2.compy/arm-diags-data/"
data_path = "/Users/zhang40/Documents/ACME/ARMDiags/data/h6_file/"
out_path = "/Users/zhang40/Documents/ACME/ARMDiags/data/processed_h6_file/"
Path(out_path).mkdir(parents=True, exist_ok=True)

# time_range = "000101_000112"
time_range = "198501_201412"
# time_range = "199601_201012"
p = Path(data_path)
# cmd = "ncrcat -h " + data_path + "*h1*nc " + data_path + "armsites_all_time.nc"
# cmd = "ncrcat -h " + data_path + "*h4*nc " + out_path + "armsites_all_time.nc"
# Below has not been tested.
cmd = (
    "ncrcat -h -d time,'1985-01-01 0:00:0.0','2014-12-31 23:59:0.0' "
    + data_path
    # + "*h4*nc "
    + "*h6*nc "
    + out_path
    + "armsites_{}.nc".format(time_range)
)
print(cmd)
os.popen(cmd).readlines()
print("create time-series")

# filename = out_path + "armsites_all_time.nc"
filename = out_path + "armsites_{}.nc".format(time_range)
print(filename)
fin = cdms2.open(filename)
variables = [
    "CLOUD",
    "CONCLD",
    "CLDICE",
    "CLDLIQ",
    "PRECT",
    "TMQ",
    "TREFHT",
    "QREFHT",
    "CLDTOT",
    "LHFLX",
    "SHFLX",
    "FLDS",
    "FSDS",
    "FLNS",
    "FSNS",
    "FLNSC",
    "FSDSC",
    "FSNSC",
    "AODVIS",
    "AODABS",
    "PS",
    "num_a1",  # Accumu mode aerosol concentration (1/kg) at lowest level
    "num_a2",  # Aitken mode aerosol concentration (1/kg) at lowest level
    "num_a3",  # Coarse mode aerosol concentration (1/kg) at lowest level
    "so4_a1",  # Accumu mode SO4 mass conc. (kg/kg) at lowest level
    "so4_a2",  # Aitken mode SO4 mass conc. (kg/kg) at lowest level
    "CCN3",  # CCN 0.1%SS concentration (1/CC) at lowest level
    "CCN4",  # CCN 0.2%SS concentration (1/CC) at lowest level
    "CCN5",  # CCN 0.5%SS concentration (1/CC) at lowest level
]
sites_info = {
    "sgpc1": [262.5, 36.6],
    "nsac1": [204.6, 71.3],
    "twpc1": [147.4, -2.0],
    "twpc2": [166.9, -0.5],
    "twpc3": [130.9, -12.4],
    "enac1": [331.97, 39.09],
}
sites = ["sgpc1", "nsac1", "twpc1", "twpc2", "twpc3", "enac1"]


for site in sites:
    if sites_info[site][1] > 0:
        lon_lat = str(sites_info[site][0]) + "e_" + str(sites_info[site][1]) + "n"
    else:
        lon_lat = str(sites_info[site][0]) + "e_" + str(abs(sites_info[site][1])) + "s"

    for variable in variables:
        fout_1 = cdms2.open(
            out_path + variable + "_" + site + "_" + time_range + ".nc", "w"
        )
        var_name = variable + "_" + lon_lat

        print(var_name)
        var = fin(var_name, squeeze=1)
        var_time = var.getTime()

        if " 0000-" in var_time.units:
            units = var_time.units
            fakeUnits = units.replace(
                " 0000-", " 0001-"
            )  # valid units for cdtime https://github.com/CDAT/cdms/issues/334
            var_time.units = fakeUnits
        var.id = variable
        lat = fin("lat_" + lon_lat)
        lat = lat[:]
        lat = cdms2.createVariable(lat)
        lat.id = "lat"
        lat.long_name = "latitude"
        lat.units = "degrees_north"
        lon = fin("lon_" + lon_lat)
        lon = lon[:]
        lon = cdms2.createVariable(lon)
        lon.id = "lon"
        lon.long_name = "longitude"
        lon.units = "degrees_east"
        fout_1.write(var)
        fout_1.write(lat)
        fout_1.write(lon)
        if var.getLevel():
            var1 = fin("PS" + "_" + lon_lat, squeeze=1)
            var1.id = "PS"
            var2 = fin("P0")
            var2.id = "P0"
            var3 = fin("hyam")
            var3.id = "hyam"
            var4 = fin("hybm")
            var4.id = "hybm"
            fout_1.write(var1)
            fout_1.write(var2)
            fout_1.write(var3)
            fout_1.write(var4)
        fout_1.close()
