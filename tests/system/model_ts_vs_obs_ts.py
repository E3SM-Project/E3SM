reference_data_path = '../testing_data/obs_for_e3sm_diags/time-series'
test_data_path = '../testing_data/test_model_data_for_e3sm_diags/time-series/E3SM_v1'

test_name = 'E3SM_v1'  # Needed so it can be shown on the plots.
results_dir = 'model_ts_vs_obs_ts_results'

ref_timeseries_input = True
ref_start_yr = 1980
ref_end_yr = 1980

test_timeseries_input = True
test_start_yr = 1850
test_end_yr = 1850

# We define the sets b/c we don't currently have any time-series data for cosp_histogram.
sets = ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar']

backend = "mpl"
debug = True
