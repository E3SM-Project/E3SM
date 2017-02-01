case_id = 'set5_ANN_T200_ECMWF'
#reference_data_path = '/space1/test_data/obs_for_diagnostics/'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'
#reference_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
#test_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
reference_data_set = 'MERRA_ANN_climo.nc'  # observation
test_data_set = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model
variables = 'T'
plev = 850
season = 'ANN'
output_file = 'test_lev.png'

# VCS options
main_title = 'this changes in the code'

test_name = '1850_alpha6_01 (yrs0070-0099)'
test_title = 'Model'
#test_colormap = ''
#test_units = 'K'

#reference_name = 'ECMWF (yrs unknown)'
reference_title = 'Observation'
#reference_colormap = ''
#reference_units = 'K'

diff_name = ''
diff_title = 'Model - Observation'
diff_colormap = 'bl_to_darkred'
#diff_units = 'K'

arrows = True
