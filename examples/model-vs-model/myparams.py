# Location of the ref data.
reference_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
# Name of the ref model data, used to find the climo files.
ref_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'
# An optional, shorter name to be used instead of the ref_name.
short_ref_name = 'Ref: beta0.F1850COSP_ne30'

# Location of the test data.
test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
# Name of the test model data, used to find the climo files.
test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
# An optional, shorter name to be used instead of the test_name.
short_test_name = 'Test: beta0_FC5COSP_ne30'

# What plotsets to run the diags on.
sets = ['lat_lon']
# Name of the folder where the results are stored.
results_dir = 'model_to_model'
# This parameter modifies the software to accommodate model vs model runs.
# The default setting for run_type is 'model_vs_obs'.
run_type = 'model_vs_model' 

# Below are more optional arguments.

# 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
backend = 'mpl'
# Title of the difference plots.
diff_title = 'Test Model - Ref Model'
# For running with multiprocessing.
multiprocessing = True
num_workers = 16
