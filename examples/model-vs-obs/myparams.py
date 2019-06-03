# Location of the data.
reference_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
# Name of the test model data, used to find the climo files.
test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
# An optional, shorter name to be used instead of the test_name.
short_test_name = 'beta0.FC5COSP.ne30'

# What plotsets to run the diags on.
sets = ['lat_lon']
# Name of the folder where the results are stored.
results_dir = 'era_tas_land'

# Below are more optional arguments.

# 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
backend = 'mpl'
# Title of the difference plots.
diff_title = 'Model - Obs.'
# Save the netcdf files for each of the ref, test, and diff plot.
save_netcdf = True
