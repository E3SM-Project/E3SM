case_id = 'set5_ANN_PRECT_TRMM'
#reference_data_path = '/space1/test_data/obs_for_diagnostics/'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'
#reference_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
#test_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
reference_data_set = 'TRMM_ANN_climo.nc'  # observation
test_data_set = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model

variables = 'PRECT'
season = 'ANN'
regrid_tool = 'esmf'
regrid_method = 'linear'

output_file = 'test'

# VCS options
main_title = 'PRECT ANN'

test_name = '1850_alpha6_01 (yrs0070-0099)'
test_title = 'Model'
test_colormap = ''
#test_levels = [0, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
test_levels = []
test_units = 'mm/day'

reference_name = 'TRMM'# (yrs1979-2009)'
reference_title = 'Observation'
reference_colormap = ''
#reference_levels = [0, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
reference_levels = []
reference_units = 'mm/day'

diff_name = ''
diff_title = 'Model - Observation'
diff_colormap = 'bl_to_darkred'
#diff_levels = [-6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6]
diff_levels = []
diff_units = 'mm/day'

canvas_size_w = 1212
canvas_size_h = 1628
