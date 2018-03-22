# Used by Jenkins on ACME1

# Example run command:
# acme_diags_driver.py -p all_sets_nightly.py --backend vcs --results_dir /var/www/acme/acme-diags/e3sm_diags_jenkins/conda_create_all_sets_vcs_`date +%Y.%m.%d-%H:%M:%S`
reference_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'

test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

# Add both of these parameters with the command line interface
# backend = ''
# results_dir = ''

multiprocessing = True
num_workers = 24

