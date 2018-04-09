reference_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
ref_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'
reference_name = 'ref: beta0.F1850COSP_ne30'

test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
short_test_name = 'test: beta0_FC5COSP_ne30'

backend = 'mpl'

diff_title = 'test mod - ref mod'
results_dir = 'model_to_model'

sets = ['lat_lon'] #without specifiying this line, all sets will run

run_type = 'model_vs_model' #This parameter modifies the look of viewer to accomodate model vs model run. Default setting for run_type is 'model_vs_obs'

#For running multiprocessing
multiprocessing = True
num_workers = 16
