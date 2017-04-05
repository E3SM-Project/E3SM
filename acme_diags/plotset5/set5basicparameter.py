#reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'

test_name = '20160520.A_WCYCL1850.ne30'

regrid_tool = 'esmf'
regrid_method = 'linear'

backend = 'vcs'

diff_title = 'test - reference'
diff_colormap = 'bl_to_darkred'

# output_format = ['png', 'pdf']

derived_variables = {
    'TEST': [
        (['test'], lambda test: 'test')
    ]
}
