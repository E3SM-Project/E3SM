import os
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from acme_diags.run import runner

param = CoreParameter()

param.reference_data_path = (
    "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/"
)
param.test_data_path = "/global/cfs/cdirs/e3sm/acme_diags/test_model_data_for_acme_diags/climatology/"
param.test_name = "20161118.beta0.FC5COSP.ne30_ne30.edison"
param.seasons = ["ANN"]

# Name of the folder where the results are stored.
# Change `prefix` to use your directory.
# prefix = '/global/cfs/cdirs/e3sm/www/<your directory>/examples'
param.results_dir = os.path.join(prefix, "ex6_zonal_mean_2d_and_lat_lon_demo")

# Uncomment the two lines below to just
# run the diags with T and PRECT.
# param.selectors += ['variables']
# param.variables = ['T', 'PRECT']

# The new changes are below.
zonal_mean_2d_param = ZonalMean2dParameter()
zonal_mean_2d_param.plevs = [10.0, 20.0, 30.0]

runner.sets_to_run = ["zonal_mean_2d", "lat_lon"]
runner.run_diags([param, zonal_mean_2d_param])
