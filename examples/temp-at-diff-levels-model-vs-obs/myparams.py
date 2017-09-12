reference_data_path = '/space1/test_data/reanalysis_data/ERA-Interim/ta/climos/'
ref_name = 'ERA-Interim_ta'

test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'
test_name = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01'

backend = 'vcs'
results_dir = 'Jerry_ra_ta'

sets = ['lat_lon']


def rename(new_name):
    '''Given the new name, just return it.'''
    return new_name


derived_variables = {
    'T': {
        ('ta'): rename
    }
}
