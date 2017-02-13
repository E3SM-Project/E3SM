reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'

#reference_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
#test_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
test_data_set = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model

regrid_tool = 'esmf'
regrid_method = 'linear'

output_file = 'test.png'

test_title = 'Model'
reference_title = 'Observation'
diff_title = 'Model - Observation'
diff_colormap = 'bl_to_darkred'

canvas_size_w = 1212
canvas_size_h = 1628
