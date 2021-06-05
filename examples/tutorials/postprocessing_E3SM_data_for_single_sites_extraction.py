import glob
import os
from subprocess import call
from typing import Dict, Tuple

import cdms2

# A script to extract time series from global E3SM output on native ne30 or ne120 output and convert to per-variable, per-site netcdf files as input for ARM diagostics.

data_path = "/global/cfs/cdirs/e3smpub/E3SM_simulations/theta.20190910.branch_noCNT.n438b.unc03.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run/"
case_name = (
    "theta.20190910.branch_noCNT.n438b.unc03.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG"
)

output_path = "/global/cfs/cdirs/e3sm/zhang40/armsites/tests/{}".format(case_name)
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Specify both time_range and file_list
time_range = "005601_005712"
file_list = glob.glob("{}{}.cam.h2.005[6-7]*.nc".format(data_path, case_name))


# ((node index-ne30, node index-ne120),[lon, lat], site name)
SiteInfo = Tuple[int, int, float, float, str]
SITE_DICT: Dict[str, SiteInfo] = {
    "sgp": (32068, 507365, 262.5, 36.6, "97.5W 36.6N Oklahoma ARM"),
    "nsa": (46132, 737479, 204.6, 71.3, "156.6W 71.3N Barrow ARM"),
    "twpc1": (20370, 321269, 147.4, -2.0, "147.4E 2.S Manus ARM"),
    "twpc2": (20435, 324740, 166.9, -0.5, "166.9E 0.5S Nauru ARM"),
    "twpc3": (11057, 172748, 130.9, -12.4, "130.9E 12.4S Darwin ARM"),
}

# Specify grid (ne30 or ne120)
res = "ne30"
if res == "ne30":
    node_ind = 0
elif res == "ne120":
    node_ind = 1

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
]

for site, SiteInfo in SITE_DICT.items():
    for file in file_list:
        file_name = os.path.basename(file)[:-3]
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
