#reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'

#reference_data_path = '/Users/shaheen2/obs_for_diagnostics/'
#test_data_path = '/Users/shaheen2/ACME_simulations/'
#test_data_set = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model
test_name = '20160520.A_WCYCL1850.ne30'


regrid_tool = 'esmf'
regrid_method = 'linear'

test_title = 'Model'
reference_title = 'Observation'
diff_title = 'Model - Observation'
diff_colormap = 'bl_to_darkred'

canvas_size_w = 1212
canvas_size_h = 1628
