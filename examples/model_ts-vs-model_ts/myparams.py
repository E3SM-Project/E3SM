# Location of the data.
test_data_path = '/p/user_pub/work/E3SM/1_0/historical_H1/1deg_atm_60-30km_ocean/atmos/129x256/time-series/mon/ens1/v1/'
reference_data_path = '/p/user_pub/work/E3SM/1_0/historical_H1/1deg_atm_60-30km_ocean/atmos/129x256/time-series/mon/ens1/v1/'

# Set this parameter to True.
# By default, e3sm_diags expects the test data to be climo data.
test_timeseries_input = True
# Years to slice the test data, base this off the years in the filenames.
test_start_yr = '2011'
test_end_yr = '2013'

# Set this parameter to True.
# By default, e3sm_diags expects the ref data to be climo data.
ref_timeseries_input = True
# Years to slice the ref data, base this off the years in the filenames.
ref_start_yr = '1850'
ref_end_yr = '1852'

# When running with time-series data, you don't need to specify the name of the data.
# But you should, otherwise nothing is displayed when the test/ref name is needed.
short_test_name = 'historical_H1'
short_ref_name = 'historical_H1'

# This parameter modifies the software to accommodate model vs model runs.
# The default setting for run_type is 'model_vs_obs'.
run_type = 'model_vs_model'
# Name of the folder where the results are stored.
results_dir = 'modTS_vs_modTS_3years'

# Below are more optional arguments.

# What plotsets to run the diags on.
# If not defined, then all available sets are used. 
sets = ['lat_lon']
# What seasons to run the diags on.
# If not defined, diags is ran on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
seasons = ['ANN']
# Title of the difference plots.
diff_title = 'Model (2011-2013) - Model (1850-1852)'
# For running with multiprocessing.
multiprocessing = True
num_workers = 24
