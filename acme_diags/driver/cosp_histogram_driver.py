from __future__ import print_function

import os
import cdms2
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean
from acme_diags.driver import utils


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    metrics_dict = {}
    metrics_dict['ref'] = {
        'min': min_cdms(ref),
        'max': max_cdms(ref),
        'mean': mean(ref)
    }
    metrics_dict['test'] = {
        'min': min_cdms(test),
        'max': max_cdms(test),
        'mean': mean(test)
    }

    metrics_dict['diff'] = {
        'min': min_cdms(diff),
        'max': max_cdms(diff),
        'mean': mean(diff)
    }
    metrics_dict['misc'] = {
        'rmse': rmse(test_regrid, ref_regrid),
        'corr': corr(test_regrid, ref_regrid)
    }

    return metrics_dict


def run_diag(parameter):
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, 'ref_name', '')
    regions = parameter.regions

    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)    

    for season in seasons:
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data, season)
        parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data, season)

        # Get land/ocean fraction for masking.
        try:
            land_frac = test_data.get_variable('LANDFRAC', season)
            ocean_frac = test_data.get_variable('OCNFRAC', season)
        except:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            with cdms2.open(mask_path) as f:
                land_frac = f('LANDFRAC')
                ocean_frac = f('OCNFRAC')

        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var

            mv1 = test_data.get_variable(var, season)
            mv2 = ref_data.get_variable(var, season)

            parameter.viewer_descr[var] = mv1.long_name if hasattr(
                mv1, 'long_name') else 'No long_name attr in test data.'

            for region in regions:
                print("Selected region: {}".format(region))

                mv1_domain, mv2_domain = utils.general.select_region(
                    region, mv1, mv2, land_frac, ocean_frac, parameter)

                parameter.output_file = '-'.join(
                    [ref_name, var, season, region])
                parameter.main_title = str(' '.join([var, season, region]))

                mv1_domain_mean = mean(mv1_domain)
                mv2_domain_mean = mean(mv2_domain)
                diff = mv1_domain_mean - mv2_domain_mean

                mv1_domain_mean.id = var
                mv2_domain_mean.id = var
                diff.id = var
                
                parameter.backend = 'mpl'  # For now, there's no vcs support for this set.
                plot(parameter.current_set, mv2_domain_mean,
                     mv1_domain_mean, diff, {}, parameter)
                utils.general.save_ncfiles(parameter.current_set,
                                   mv1_domain, mv2_domain, diff, parameter)

    return parameter
