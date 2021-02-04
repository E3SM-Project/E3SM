import os

from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.run import runner

param = CoreParameter()

# Location of the ref data.
param.reference_data_path = "/global/cfs/cdirs/e3sm/acme_diags/test_model_data_for_acme_diags/climatology/"
# Name of the ref model data, used to find the climo files.
param.ref_name = "20161118.beta0.F1850COSP.ne30_ne30.edison"
# An optional, shorter name to be used instead of the ref_name.
param.short_ref_name = "Ref: beta0.F1850COSP_ne30"

# Location of the test data.
param.test_data_path = "/global/cfs/cdirs/e3sm/acme_diags/test_model_data_for_acme_diags/climatology"
# Name of the test model data, used to find the climo files.
param.test_name = "20161118.beta0.FC5COSP.ne30_ne30.edison"
# An optional, shorter name to be used instead of the test_name.
param.short_test_name = "Test: beta0_FC5COSP_ne30"

# What plotsets to run the diags on.
param.sets = ["lat_lon"]
# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
prefix = "/global/cfs/cdirs/e3sm/www/<your directory>/examples"
param.results_dir = os.path.join(prefix, "ex4_model_to_model")
# This parameter modifies the software to accommodate model vs model runs.
# The default setting for run_type is 'model_vs_obs'.
param.run_type = "model_vs_model"

# Below are more optional arguments.

# 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
param.backend = "mpl"
# Title of the difference plots.
param.diff_title = "Test Model - Ref Model"
# For running with multiprocessing.
# multiprocessing = True
# num_workers = 32

runner.sets_to_run = ["lat_lon"]
runner.run_diags([param])
