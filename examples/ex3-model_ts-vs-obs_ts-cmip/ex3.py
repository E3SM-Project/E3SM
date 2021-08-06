import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# Location of the data.
param.test_data_path = "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/CESM1-CAM5_cmip/"
param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/time-series/"
)

# Set this parameter to True.
# By default, e3sm_diags expects the test data to be climo data.
param.test_timeseries_input = True
# Years to slice the test data, base this off the years in the filenames.
param.test_start_yr = "2003"
param.test_end_yr = "2004"

# Set this parameter to True.
# By default, e3sm_diags expects the ref data to be climo data.
param.ref_timeseries_input = True
# Years to slice the ref data, base this off the years in the filenames.
param.ref_start_yr = "2003"
param.ref_end_yr = "2004"

# When running with time-series data, you don't need to specify the name of the data.
# But you should, otherwise nothing is displayed when the test/ref name is needed.
param.short_test_name = "CESM1-CAM5-historical"
# param.short_ref_name = ''

# This parameter modifies the software to accommodate model vs model runs.
# The default setting for run_type is 'model_vs_obs'.
# param.run_type = 'model_vs_model'
# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/<your directory>/examples"
param.results_dir = os.path.join(prefix, "ex3_modTS_CMIP_vs_obs_2years")

# Below are more optional arguments.

# What plotsets to run the diags on.
# If not defined, then all available sets are used.
param.sets = ["lat_lon"]
# What seasons to run the diags on.
# If not defined, diags are ran on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
param.seasons = ["ANN"]
# Title of the difference plots.
param.diff_title = "model(2003-2004 yr_avg) - obs(2003-2004 yr_avg)"
# For running with multiprocessing.
# param.multiprocessing = True
# param.num_workers = 32

runner.sets_to_run = ["lat_lon"]
runner.run_diags([param])
