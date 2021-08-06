import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# Location of the data.
param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/"
)
param.test_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/climatology/"
)
# Name of the test model data, used to find the climo files.
param.test_name = "20161118.beta0.FC5COSP.ne30_ne30.edison"
# An optional, shorter name to be used instead of the test_name.
param.short_test_name = "beta0.FC5COSP.ne30"

# What plotsets to run the diags on.
param.sets = ["lat_lon"]
# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/<your directory>/examples"
param.results_dir = os.path.join(prefix, "ex5_model_to_obs")

# Below are more optional arguments.

# 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
param.backend = "mpl"
# Title of the difference plots.
param.diff_title = "Model - Obs."
# Save the netcdf files for each of the ref, test, and diff plot.
param.save_netcdf = True
# For running with multiprocessing.
# param.multiprocessing = True
# param.num_workers = 32

runner.sets_to_run = ["lat_lon"]
runner.run_diags([param])
