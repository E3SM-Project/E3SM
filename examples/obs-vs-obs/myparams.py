# Location of the ref data.
reference_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
# Name of the ref obs data, used to find the climo files.
ref_name = 'ceres_ebaf_toa_v2.8'

# Location of the test data.
test_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
# Name of the test obs data, used to find the climo files.
test_name = 'ceres_ebaf_toa_v4.0'

# Name of the folder where the results are stored.
results_dir = 'obs_vs_obs'
# What plotsets to run the diags on.
sets = ['lat_lon']

# Below are more optional arguments.

# What seasons to run the diags on.
# If not defined, diags is ran on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
seasons = ['ANN']
# 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
backend = 'mpl'
