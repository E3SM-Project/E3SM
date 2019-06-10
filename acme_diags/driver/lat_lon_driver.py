from __future__ import print_function

import os
import json
import cdms2
import MV2
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean, std
from acme_diags.driver import utils


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    metrics_dict = {}
    metrics_dict['ref'] = {
        'min': float(min_cdms(ref)),
        'max': float(max_cdms(ref)),
        'mean': float(mean(ref))
    }
    metrics_dict['ref_regrid'] = {
        'min': float(min_cdms(ref_regrid)),
        'max': float(max_cdms(ref_regrid)),
        'mean': float(mean(ref_regrid)),
        'std': float(std(ref_regrid))
    }
    metrics_dict['test'] = {
        'min': float(min_cdms(test)),
        'max': float(max_cdms(test)),
        'mean': float(mean(test))
    }
    metrics_dict['test_regrid'] = {
        'min': float(min_cdms(test_regrid)),
        'max': float(max_cdms(test_regrid)),
        'mean': float(mean(test_regrid)),
        'std': float(std(test_regrid))
    }
    metrics_dict['diff'] = {
        'min': float(min_cdms(diff)),
        'max': float(max_cdms(diff)),
        'mean': float(mean(diff))
    }
    metrics_dict['misc'] = {
        'rmse': float(rmse(test_regrid, ref_regrid)),
        'corr': float(corr(test_regrid, ref_regrid))
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

            # Special case, cdms didn't properly convert mask with fill value
            # -999.0, filed issue with Denis.
            if ref_name == 'WARREN':
                # This is cdms2 fix for bad mask, Denis' fix should fix this.
                mv2 = MV2.masked_where(mv2 == -0.9, mv2)
            # The following should be moved to a derived variable.
            if ref_name == 'AIRS':
                # This is cdms2 fix for bad mask, Denis' fix should fix this.
                mv2 = MV2.masked_where(mv2 > 1e+20, mv2)
            if ref_name == 'WILLMOTT' or ref_name == 'CLOUDSAT':
                # This is cdms2 fix for bad mask, Denis' fix should fix this.
                mv2 = MV2.masked_where(mv2 == -999., mv2)

                # The following should be moved to a derived variable.
                if var == 'PRECT_LAND':
                    days_season = {'ANN': 365, 'DJF': 90,
                                   'MAM': 92, 'JJA': 92, 'SON': 91}
                    # mv1 = mv1 * days_season[season] * 0.1 # following AMWG
                    # Approximate way to convert to seasonal cumulative
                    # precipitation, need to have solution in derived variable,
                    # unit convert from mm/day to cm.
                    mv2 = mv2 / days_season[season] / \
                        0.1  # Convert cm to mm/day instead.
                    mv2.units = 'mm/day'

            # For variables with a z-axis.
            if mv1.getLevel() and mv2.getLevel():
                plev = parameter.plevs
                print('Selected pressure level: {}'.format(plev))

                mv1_p = utils.general.convert_to_pressure_levels(mv1, plev, test_data, var, season)
                mv2_p = utils.general.convert_to_pressure_levels(mv2, plev, ref_data, var, season)

                # Select plev.
                for ilev in range(len(plev)):
                    mv1 = mv1_p[ilev, ]
                    mv2 = mv2_p[ilev, ]

                    for region in regions:
                        print("Selected region: {}".format(region))

                        mv1_domain, mv2_domain = utils.general.select_region(
                            region, mv1, mv2, land_frac, ocean_frac, parameter)

                        parameter.output_file = '-'.join(
                            [ref_name, var, str(int(plev[ilev])), season, region])
                        parameter.main_title = str(
                            ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))

                        # Regrid towards the lower resolution of the two
                        # variables for calculating the difference.
                        mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
                            mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                        diff = mv1_reg - mv2_reg
                        metrics_dict = create_metrics(
                            mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)

                        # Saving the metrics as a json.
                        metrics_dict['unit'] = mv1_reg.units
                        fnm = os.path.join(utils.general.get_output_dir(
                            parameter.current_set, parameter), parameter.output_file + '.json')
                        with open(fnm, 'w') as outfile:
                             json.dump(metrics_dict, outfile)
                        # Get the filename that the user has passed in and display that.
                        # When running in a container, the paths are modified.
                        fnm = os.path.join(utils.general.get_output_dir(parameter.current_set,
                            parameter, ignore_container=True), parameter.output_file + '.json')
                        print('Metrics saved in: ' + fnm)

                        parameter.var_region = region
                        plot(parameter.current_set, mv2_domain,
                             mv1_domain, diff, metrics_dict, parameter)
                        utils.general.save_ncfiles(
                            parameter.current_set, mv1_domain, mv2_domain, diff, parameter)


            # For variables without a z-axis.
            elif mv1.getLevel() is None and mv2.getLevel() is None:
                for region in regions:
                    print("Selected region: {}".format(region))

                    mv1_domain, mv2_domain = utils.general.select_region(
                        region, mv1, mv2, land_frac, ocean_frac, parameter)

                    parameter.output_file = '-'.join(
                        [ref_name, var, season, region])
                    parameter.main_title = str(' '.join([var, season, region]))

                    # Regrid towards the lower resolution of the two
                    # variables for calculating the difference.
                    mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
                        mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                    # Special case.
                    if var == 'TREFHT_LAND' or var == 'SST':
                        if ref_name == 'WILLMOTT':
                            mv2_reg = MV2.masked_where(
                                mv2_reg == mv2_reg.fill_value, mv2_reg)

                        land_mask = MV2.logical_or(mv1_reg.mask, mv2_reg.mask)
                        mv1_reg = MV2.masked_where(land_mask, mv1_reg)
                        mv2_reg = MV2.masked_where(land_mask, mv2_reg)

                    diff = mv1_reg - mv2_reg
                    metrics_dict = create_metrics(
                        mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)

                    # Saving the metrics as a json.
                    metrics_dict['unit'] = mv1_reg.units
                    fnm = os.path.join(utils.general.get_output_dir(
                        parameter.current_set, parameter), parameter.output_file + '.json')
                    with open(fnm, 'w') as outfile:
                            json.dump(metrics_dict, outfile)
                    # Get the filename that the user has passed in and display that.
                    # When running in a container, the paths are modified.
                    fnm = os.path.join(utils.general.get_output_dir(parameter.current_set,
                        parameter, ignore_container=True), parameter.output_file + '.json')
                    print('Metrics saved in: ' + fnm)

                    parameter.var_region = region
                    plot(parameter.current_set, mv2_domain,
                         mv1_domain, diff, metrics_dict, parameter)
                    utils.general.save_ncfiles(parameter.current_set,
                                       mv1_domain, mv2_domain, diff, parameter)

            else:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting.")

    return parameter
