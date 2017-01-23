case_id = 'set5_ANN_PRECT_GPCP'
#reference_data_path = '/space1/test_data/obs_for_diagnostics/'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'
reference_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
test_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
reference_data_set = 'GPCP_ANN_climo.nc'  # observation
test_data_set = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model
variables = 'PRECT'
season = 'ANN'
output_file = 'test'

# All of the metrics that Chris wants
# ONLY HAVE THE 'MOST IMPORTANT' ONES
reference_name = 'GPCP_1979-2009' # appears on the top right of the graph
test_name = 'GPCP_1979-2009'
diff_name = 'GPCP_1979-2009'

main_title = 'PRECT ANN'
reference_title = 'observation'
test_title = 'reference'
diff_title = 'model - observation'

canvas_size_w = 1212
canvas_size_h = 1628
output_file_format = 'png'
logo = False
reference_levels = []
reference_colormap = ''

test_levels = []
test_colormap = ''

diff_levels = []
diff_colormap = ''

regrid_method = ''
regrid_tool = ''
