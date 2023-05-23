import glob
import os
from subprocess import call
from typing import Dict, Tuple

import cdms2

# A script to extract time series from global E3SM output on native ne30 or ne120 output and convert to per-variable, per-site netcdf files as input for ARM diagostics.

# case_name = "theta.20190910.branch_noCNT.n438b.unc03.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG"
# data_path = "/global/cfs/cdirs/e3smpub/E3SM_simulations/{}/run/".format(case_name)
#
# Specify both time_range and file_list
# time_range = "005601_005712"
# file_list = glob.glob("{}{}.cam.h4.005[6-7]*.nc".format(data_path, case_name))  #h4 tape outputs 3hourly data
#
#
# Specify grid (ne30 or ne120)
# res = "ne120"
# pg2 = False

# E3SM v2 output with ne30gp2
case_name = "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
data_path = "/global/cfs/cdirs/e3smdata/zppy_complete_run_nersc_output/{}/archive/atm/hist/".format(
    case_name
)
time_range = "005101_005212"
file_list = glob.glob(
    "{}{}.eam.h4.005[1-2]*.nc".format(data_path, case_name)
)  # h4 tape outputs 3hourly data

# Specify grid (ne30 or ne120)
res = "ne30"
pg2 = True


if res == "ne30":
    node_ind = 0
elif res == "ne120":
    node_ind = 1
# ((node index-ne30, node index-ne120),[lon, lat], site name)
SiteInfo = Tuple[int, int, float, float, str]
if pg2:
    SITE_DICT: Dict[str, SiteInfo] = {
        "sgpc1": (14089, 225317, 262.5, 36.6, "97.5W 36.6N Oklahoma ARM"),
        "nsac1": (20450, 327558, 204.6, 71.3, "156.6W 71.3N Barrow ARM"),
        "twpc1": (8896, 142147, 147.4, -2.0, "147.4E 2.S Manus ARM"),
        "twpc2": (8923, 143689, 166.9, -0.5, "166.9E 0.5S Nauru ARM"),
        "twpc3": (4795, 76297, 130.9, -12.4, "130.9E 12.4S Darwin ARM"),
    }
else:
    SITE_DICT: Dict[str, SiteInfo] = {
        "sgpc1": (32068, 507365, 262.5, 36.6, "97.5W 36.6N Oklahoma ARM"),
        "nsac1": (46132, 737479, 204.6, 71.3, "156.6W 71.3N Barrow ARM"),
        "twpc1": (20370, 321269, 147.4, -2.0, "147.4E 2.S Manus ARM"),
        "twpc2": (20435, 324740, 166.9, -0.5, "166.9E 0.5S Nauru ARM"),
        "twpc3": (11057, 172748, 130.9, -12.4, "130.9E 12.4S Darwin ARM"),
    }


# Node index for pg2 created by Wuyin Lin with an IDL script. Will incorporate a python script to locate the nearest columns for any model grids.
# Grid: ne30pg2
# -------------
# SGP  lon      262.500  lat      36.6000  nearest columns (0-based)
#            263.24867       36.559257       14089
#
# NSA  lon      204.600  lat      71.3000  nearest columns (0-based)
#
#            205.03112       71.083494       20450
# TWPC1  lon      147.400  lat     -2.00000  nearest columns (0-based)
#            147.74997      -1.9033821        8896
#
# TWPC2  lon      166.900  lat    -0.500000  nearest columns (0-based)
#            167.25000     -0.73153949        8923
#
# TWPC3  lon      130.900  lat     -12.4000  nearest columns (0-based)
#            131.25111      -11.974293        4795
#
# Grid: ne120pg2
# --------------
# SGP  lon      262.500  lat      36.6000  nearest columns (0-based)
#            262.31241       36.689360      225317
#
# NSA  lon      204.600  lat      71.3000  nearest columns (0-based)
#            204.77460       71.322908      327558
#
# TWPC1  lon      147.400  lat     -2.00000  nearest columns (0-based)
#            147.56250      -2.0575710      142147
#
# TWPC2  lon      166.900  lat    -0.500000  nearest columns (0-based)
#            167.06250     -0.54822243      143689
#
# TWPC3  lon      130.900  lat     -12.4000  nearest columns (0-based)
#            131.06257      -12.443734       76297


output_path = "/global/cfs/cdirs/e3sm/zhang40/armsites/tests/{}".format(case_name)
if not os.path.exists(output_path):
    os.makedirs(output_path)

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
    "QFLX",
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

for site, SiteInfo in SITE_DICT.items():
    for file in file_list:
        file_name = os.path.basename(file)[:-3]
        print(file_name)
        print(file_name, res, SiteInfo[node_ind])
        cmd = "ncks -d ncol,{},{} {} {}/{}_{}.nc".format(
            SiteInfo[node_ind], SiteInfo[node_ind], file, output_path, file_name, site
        )
        call(cmd, shell=True)

    output_file_name = "{}_{}_all_time.nc".format(case_name, site)
    cmd = "ncrcat -h {}/{}*{}.nc {}/{}".format(
        output_path, case_name, site, output_path, output_file_name
    )
    print(cmd)
    call(cmd, shell=True)

    fin = cdms2.open("{}/{}".format(output_path, output_file_name))
    for variable in variables:
        if variable in fin.listvariable():
            print(variable)
            var = fin(variable, squeeze=1)
            fout = cdms2.open(
                "{}/{}_{}_{}.nc".format(output_path, variable, site, time_range), "w"
            )
            lat = SiteInfo[3]
            lat = cdms2.createVariable(lat)
            lat.id = "lat"
            lat.long_name = "latitude"
            lat.units = "degrees_north"
            lon = SiteInfo[2]
            lon = cdms2.createVariable(lon)
            lon.id = "lon"
            lon.long_name = "longitude"
            lon.units = "degrees_east"

            fout.write(var)
            fout.write(lat)
            fout.write(lon)
            fout.close()
    fin.close()
