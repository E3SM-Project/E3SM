import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

# Location of the ref data.
param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/"
)
# Name of the ref obs data, used to find the climo files.
param.ref_name = "ceres_ebaf_toa_v2.8"

# Location of the test data.
param.test_data_path = (
    "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/climatology/"
)
# Name of the test obs data, used to find the climo files.
param.test_name = "ceres_ebaf_toa_v4.0"

# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/<your directory>/examples"
param.results_dir = os.path.join(prefix, "ex7_obs_vs_obs")
# What plotsets to run the diags on.
param.sets = ["lat_lon"]

# Below are more optional arguments.

# What seasons to run the diags on.
# If not defined, diags is ran on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
param.seasons = ["ANN"]
# 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
param.backend = "mpl"

runner.sets_to_run = ["lat_lon"]
runner.run_diags([param])
