# Used by Jenkins on ACME1

# Example run command:
# acme_diags -p all_sets_nightly_model_vs_model.py --backend vcs --results_dir /var/www/acme/acme-diags/e3sm_diags_jenkins/conda_create_all_sets_vcs_`date +%Y.%m.%d-%H:%M:%S`
reference_data_path = "/p/cscratch/acme/data/test_model_data_for_acme_diags/"
ref_name = "20161118.beta0.F1850COSP.ne30_ne30.edison"
reference_name = "ref: beta0.F1850COSP_ne30"

test_data_path = "/p/cscratch/acme/data/test_model_data_for_acme_diags/"
test_name = "20161118.beta0.FC5COSP.ne30_ne30.edison"
short_test_name = "test: beta0_FC5COSP_ne30"

# Add both of these parameters with the command line interface
# backend = ''
# results_dir = ''

run_type = "model_vs_model"

multiprocessing = True
num_workers = 24
