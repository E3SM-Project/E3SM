import os
from acme_diags.run import runner
from acme_diags.parameter.core_parameter import CoreParameter

param = CoreParameter()

param.reference_data_path = '/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology'
param.test_data_path = '/global/cfs/cdirs/e3sm/acme_diags/test_model_data_for_acme_diags/climatology/'
param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
param.seasons = ["ANN", "JJA"]    #Default setting: seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
param.multiprocessing = True
param.num_workers = 32

#Additional parameters:
#param.short_test_name = 'beta0.FC5COSP.ne30'
#param.run_type = 'model_vs_model'
#param.diff_title = 'Difference'
#param.output_format = ['png']     
#param.output_format_subplot = ['pdf']
#param.save_netcdf = True

prefix = '/global/cfs/cdirs/e3sm/www/zhang40/tutorial2020'
param.results_dir = os.path.join(prefix, 'climo_sets_haswell')

runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
runner.run_diags([param])

