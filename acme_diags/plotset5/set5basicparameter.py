#reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
#reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
#test_data_path = '/Users/zhang40/Documents/ACME_simulations/'

reference_data_path = '/Users/shaheen2/obs_for_diagnostics/'
test_data_path = '/Users/shaheen2/ACME_simulations/'

#reference_data_path = '/home/golaz1/ACME/diagnostics/data/obs_for_diagnostics/'
#test_data_path = '/home/golaz1/ACME/diagnostics/data/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

#test_data_set = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model
test_name = '20160520.A_WCYCL1850.ne30'

regrid_tool = 'esmf'
regrid_method = 'linear'

#test_title = 'Model'
#reference_title = 'Observation'

backend = 'vcs'
diff_title = 'test - reference'
diff_colormap = 'bl_to_darkred'

output_format = ['png', 'pdf']

def my_fcn(ref, test, diff, metrics_dict, parameter):
    print(':)')

# plot = my_fcn