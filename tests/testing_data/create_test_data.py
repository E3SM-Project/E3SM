"""
Script used to create the test data in this directory.
"""

class FileAndVars:
    def __init__(self, remote_uri, local_uri, vars_to_keep, time_slice=()):
        self.remote_uri = remote_uri
        self.local_uri = local_uri
        self.vars_to_keep = vars_to_keep
        self.time_slice = time_slice

def run_cmd(cmd):
    print('+ {}'.format(cmd))
    cmd = cmd.split()
    p = subprocess.Popen(cmd)
    p.communicate()

import os
import subprocess
import cdms2

# The data is on Cooley @ ALCF.
prefix = '/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/'

print('WARNING: Downloading about 35GB of data.')
print('You might want to do this piecewise.')

files_and_vars = []

# Obs climatology files.
path = 'obs_for_e3sm_diags/climatology/GPCP_v2.3/GPCP_v2.3_ANN_climo.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path,
    ['PRECT']
)
files_and_vars.append(f_and_v)

path = 'obs_for_e3sm_diags/climatology/MERRA2/MERRA2_ANN_198001_201612_climo.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path,
    ['ta']
)
files_and_vars.append(f_and_v)

path = 'obs_for_e3sm_diags/climatology/MISRCOSP/MISRCOSP_ANN_climo.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path,
    ['CLMISR']
)
files_and_vars.append(f_and_v)

# Obs timeseries files.
path = 'obs_for_e3sm_diags/time-series/GPCP_v2.3/PRECT_197901_201712.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path.replace('2017', '1983'),  # 1979 to 1983.
    ['PRECT'],
    time_slice=(0, 5*12)  # 5 years.
)
files_and_vars.append(f_and_v)

path = 'obs_for_e3sm_diags/time-series/MERRA2/ta_198001_201612.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path.replace('2016', '1980'),  # 01-1980 to 12-1980.
    ['ta'],
    time_slice=(0, 12)  # 1 year.
)
files_and_vars.append(f_and_v)

# There's no obs timeseries data for cosp variables.

# Model climatology files.
path = 'test_model_data_for_e3sm_diags/climatology/20161118.beta0.F1850COSP.ne30_ne30.edison_ANN_climo.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path,
    ['PRECC', 'PRECL', 'T', 'CLD_MISR'],
)
files_and_vars.append(f_and_v)

# Model timeseries files.
path = 'test_model_data_for_e3sm_diags/time-series/E3SM_v1/PRECC_185001_201312.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path.replace('2013', '1854'),  # 1850 to 1854.
    ['PRECC'],
    time_slice=(0, 5*12)  # 5 years.
)
files_and_vars.append(f_and_v)

path = 'test_model_data_for_e3sm_diags/time-series/E3SM_v1/PRECL_185001_201312.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path.replace('2013', '1854'),  # 1850 to 1854.
    ['PRECL'],
    time_slice=(0, 5*12)  # 5 years.
)
files_and_vars.append(f_and_v)

path = 'test_model_data_for_e3sm_diags/time-series/E3SM_v1/T_185001_201312.nc'
f_and_v = FileAndVars(
    os.path.join(prefix, path),
    path.replace('2013', '1850'),  # 01-1850 to 12-1850.
    ['T'],
    time_slice=(0, 12)  # 1 year.
)
files_and_vars.append(f_and_v)

# There's no model timeseries data for cosp variables.



# scp all of these files to the TEMP_DIR.
TEMP_DIR = 'temp'
if not os.path.exists(TEMP_DIR):
    os.makedirs(TEMP_DIR)
"""
for f_and_v in files_and_vars:
    cmd = 'scp zshaheen@cooley.alcf.anl.gov:{} {}'.format(f_and_v.remote_uri, TEMP_DIR)
    run_cmd(cmd)
"""

for f_and_v in files_and_vars:
    # Make the local directories.
    path, file_name = '/'.join(f_and_v.local_uri.split('/')[:-1]), f_and_v.remote_uri.split('/')[-1]
    if not os.path.exists(path):
        os.makedirs(path)

    # Open the file, get the variables we need, and slice if needed.
    orig_path = os.path.join(TEMP_DIR, file_name)
    with cdms2.open(orig_path) as orig_file:
        if os.path.exists(f_and_v.local_uri):
            os.remove(f_and_v.local_uri)

        with cdms2.open(f_and_v.local_uri, 'w') as new_file:
            # Copy over the attributes.
            new_file.attributes = orig_file.attributes

            # Lat, lon, and time are considered dimensions in cdms2.
            # Those are automatically saved.
            other_vars = ['bounds_time', 'bounds_lat',
                            'bounds_lon', 'hyam', 'hybm', 'PS']
            for v in f_and_v.vars_to_keep + other_vars:
                if v not in orig_file.variables:
                    msg = 'ERROR: {} not in {}'.format(v, orig_path)
                    print(msg)
                    continue
                if f_and_v.time_slice:
                    var = orig_file(v, time=slice(f_and_v.time_slice[0], f_and_v.time_slice[1]))
                else:
                    var = orig_file(v)
                new_file.write(var)
            print('Saved {}'.format(f_and_v.local_uri))
    
# os.remove(TEMP_DIR)
