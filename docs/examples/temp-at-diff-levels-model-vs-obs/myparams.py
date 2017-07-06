reference_data_path = '/Users/shaheen2/obs_for_diagnostics/'
ref_name = 'ERA-Interim_ta'

test_data_path = '/Users/shaheen2/ACME_simulations/'
test_name = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01'

backend = 'cartopy'
diff_colormap = 'bl_to_darkred'
results_dir = 'Jerry_ra_ta'

sets = [5]

def rename(new_name):
    '''Given the new name, just return it.'''
    return new_name

derived_variables = {
    'T': {
        ('ta'): rename
    }
}

