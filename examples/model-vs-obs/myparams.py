reference_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'  #String to search for test model
short_test_name = 'beta0.FC5COSP.ne30'   # user specified name shown on plots

backend = 'mpl'                          # select 'matplotlib as backend
diff_title = 'model - obs.'
sets=['lat_lon']
results_dir = 'era_tas_land'

save_netcdf = True   #save netcdf file being plotted
