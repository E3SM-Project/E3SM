#!/usr/bin/env python3


import os
import sys

# add eamxx scripts dir to get stuff from there
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../../scripts'))
# get the pylib stuff and ensure we have xarray
from utils import _ensure_pylib_impl
_ensure_pylib_impl('xarray')
_ensure_pylib_impl('numpy')


# get the input arg which is the case dir
case_dir = sys.argv[1]
# parse the case name, which is the string after the last "/"
case_name = case_dir.split('/')[-1]

import xarray as xr

# before the restart
ds_avg = xr.open_dataset(case_dir + '/run/' + case_name + '.scream.average_5_steps.h.AVERAGE.nsteps_x5.0001-01-01-00000.nc')
ds_ins = xr.open_dataset(case_dir + '/run/' + case_name + '.scream.instant_every_step.h.INSTANT.nsteps_x1.0001-01-01-03600.nc')

import numpy as np
assert np.allclose(ds_avg.isccp_cldtot.values, ds_ins.isccp_cldtot.mean('time', skipna=True).values, equal_nan=True)

print("PASS")
sys.exit(0)
