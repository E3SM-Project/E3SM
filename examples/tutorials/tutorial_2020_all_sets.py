import os

from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.run import runner

# See https://e3sm-project.github.io/e3sm_diags/docs/html/quickguides/quick-guide-general.html
# Compy
# data_prefix = '/compyfs/e3sm_diags_data'
# html_prefix = '/compyfs/www/<username>'  # Change <username>
# Cori
data_prefix = "/global/cfs/cdirs/e3sm/e3sm_diags"
html_prefix = "/global/cfs/cdirs/e3sm/www/<username>"  # Change <username>

param = CoreParameter()

param.reference_data_path = os.path.join(data_prefix, "obs_for_e3sm_diags/climatology")
param.test_data_path = os.path.join(
    data_prefix, "test_model_data_for_acme_diags/climatology/"
)
param.test_name = "20161118.beta0.FC5COSP.ne30_ne30.edison"
param.seasons = [
    "ANN",
    "JJA",
]  # Default setting: seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

param.results_dir = os.path.join(html_prefix, "tutorial_2020_all_sets")
param.multiprocessing = True
param.num_workers = 30

# Additional parameters:
# param.short_test_name = 'beta0.FC5COSP.ne30'
# param.run_type = 'model_vs_model'
# param.diff_title = 'Difference'
# param.output_format = ['png']
# param.output_format_subplot = ['pdf']
# param.save_netcdf = True


# Set specific parameters for new sets
enso_param = EnsoDiagsParameter()
enso_param.reference_data_path = os.path.join(
    data_prefix, "obs_for_e3sm_diags/time-series/"
)
enso_param.test_data_path = os.path.join(
    data_prefix, "test_model_data_for_acme_diags/time-series/E3SM_v1/"
)
enso_param.test_name = "e3sm_v1"
enso_param.start_yr = "1990"
enso_param.end_yr = "1999"

qbo_param = QboParameter()
qbo_param.reference_data_path = os.path.join(
    data_prefix, "obs_for_e3sm_diags/time-series/"
)
qbo_param.test_data_path = os.path.join(
    data_prefix, "test_model_data_for_acme_diags/time-series/E3SM_v1/"
)
qbo_param.test_name = "e3sm_v1"
qbo_param.start_yr = "1990"
qbo_param.end_yr = "1999"

ts_param = AreaMeanTimeSeriesParameter()
ts_param.reference_data_path = os.path.join(
    data_prefix, "obs_for_e3sm_diags/time-series/"
)
ts_param.test_data_path = os.path.join(
    data_prefix, "test_model_data_for_acme_diags/time-series/E3SM_v1/"
)
ts_param.test_name = "e3sm_v1"
ts_param.start_yr = "1990"
ts_param.end_yr = "1999"

runner.sets_to_run = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "polar",
    "cosp_histogram",
    "meridional_mean_2d",
    "enso_diags",
    "qbo",
    "area_mean_time_series",
]
runner.run_diags([param, enso_param, qbo_param, ts_param])
