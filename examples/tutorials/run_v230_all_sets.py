import os

from acme_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from acme_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from acme_diags.parameter.qbo_parameter import QboParameter
from acme_diags.parameter.streamflow_parameter import StreamflowParameter
from acme_diags.run import runner

# Define data paths for obs
input_prefix = "/global/cfs/cdirs/e3sm/acme_diags"
obs_climo = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology"
obs_ts = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/time-series"

# Define data paths for test model
test_prefix = "/global/cfs/cdirs/e3sm/zhang40/postprocessing_for_e3sm_diags"
case = "20180215.DECKv1b_H1.ne30_oEC.edison"
casename = case.split(".")[1] + "." + case.split(".")[2]

climo_path = os.path.join(test_prefix, "climo/" + case + "/1980-2014/rgr")
ts_path = os.path.join(test_prefix, "monthly_ts/" + case + "/1980-2014/rgr")
dc_climo_path = os.path.join(test_prefix, "diurnal_climo/" + case + "/1980-2014/rgr")

# Define parameters for core sets: 'lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d'
param = CoreParameter()
param.reference_data_path = obs_climo
param.test_data_path = climo_path
param.test_name = case
param.test_short_name = casename

# Define results path.
prefix = "/global/cfs/cdirs/e3sm/www/zhang40/tutorials/"
param.results_dir = os.path.join(prefix, "run_v230_allsets")

param.multiprocessing = True
param.num_workers = 30
param.seasons = ["ANN", "JJA"]

# Additional parameters:
# Short version of test model name being printed as plot titles
# param.short_test_name = 'beta0.FC5COSP.ne30'
# Specify run_type. Defualt is 'model_vs_obs'
# param.run_type = 'model_vs_model'
# Specify title of the 3rd panel plot. Defualt is 'Model - Observation'
# param.diff_title = 'Difference'
# Save subplots as pdf files.
# param.output_format_subplot = ['pdf']
# Save netcdf files being plotted.
# param.save_netcdf = True


# Set specific parameters for new sets
qbo_param = QboParameter()
qbo_param.reference_data_path = obs_ts
qbo_param.test_data_path = ts_path
qbo_param.test_name = casename
qbo_param.start_yr = "1980"
qbo_param.end_yr = "2014"


dc_param = DiurnalCycleParameter()
dc_param.reference_data_path = obs_climo
dc_param.test_data_path = dc_climo_path
dc_param.test_name = case
dc_param.short_test_name = casename
# Plotting diurnal cycle amplitude on different scales. Default is True
dc_param.normalize_test_amp = False

enso_param = EnsoDiagsParameter()
enso_param.reference_data_path = obs_ts
enso_param.test_data_path = ts_path
enso_param.test_name = casename
enso_param.start_yr = "1980"
enso_param.end_yr = "2014"

ts_param = AreaMeanTimeSeriesParameter()
ts_param.reference_data_path = obs_ts
ts_param.test_data_path = ts_path
ts_param.test_name = casename
ts_param.start_yr = "1980"
ts_param.end_yr = "2014"

streamflow_param = StreamflowParameter()
streamflow_param.reference_data_path = obs_ts
streamflow_param.test_data_path = ts_path
streamflow_param.test_start_yr = "1980"
streamflow_param.test_end_yr = "2014"
# Streamflow gauge station data range from year 1986 to 1995
streamflow_param.ref_start_yr = "1986"
streamflow_param.ref_end_yr = "1995"
#
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
    "diurnal_cycle",
    "streamflow",
]
runner.run_diags([param, enso_param, qbo_param, ts_param, dc_param, streamflow_param])
