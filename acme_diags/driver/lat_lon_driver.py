from __future__ import print_function

import os
import sys
import json
import traceback
import cdms2
import MV2
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean, std
from acme_diags.driver import utils
from acme_diags.driver.utils import get_output_dir


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """
    Creates the mean, max, min, rmse, corr in a dictionary
    """
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

def _get_land_or_ocean_frac(test_file, frac_type):
    """
    Get the land or ocean fraction for masking.
    """
    try:
        frac = test_file(frac_type)
    except:
        mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
        with cdms2.open(mask_path) as f:
            frac = f(frac_type)

    return frac

def _handle_special_cases_for_ref(var, ref_name, season, ref):
    """
    Handle the special cases for the reference data
    """
    if ref_name == 'WARREN':
        # This is a cdms fix for a bad mask.
        ref = MV2.masked_where(ref == -0.9, ref)

    elif ref_name == 'AIRS':
        # This is a cdms fix for a bad mask.
        ref = MV2.masked_where(ref > 1e+20, ref)

    # cdms didn't properly convert mask with fill value -999.0
    elif ref_name == 'WILLMOTT' or ref_name == 'CLOUDSAT':
        # This is a cdms fix for a bad mask.
        ref = MV2.masked_where(ref == -999.0, ref)

        # following should move to derived variable
        if var == 'PRECT_LAND':
            days_season = {'ANN': 365, 'DJF': 90,
                            'MAM': 92, 'JJA': 92, 'SON': 91}
            # approximate way to convert to seasonal cumulative
            # precipitation, need to have solution in derived variable,
            # unit convert from mm/day to cm
            ref = ref / days_season[season] / \
                0.1  # convert cm to mm/day instead
            ref.units = 'mm/day'
            
    return ref

def _convert_to_pressure_levels(mv, mv_file, plevs):
    """
    Given either test or reference data with a z-axis,
    convert to the desired pressure levels.
    """
    mv_plv = mv.getLevel()
    # var(time,lev,lon,lat) convert from hybrid level to pressure
    if mv_plv.long_name.lower().find('hybrid') != -1:
        hyam = mv_file('hyam')
        hybm = mv_file('hybm')
        ps = mv_file('PS')  # Pa
        mv_p = utils.hybrid_to_plevs(mv, hyam, hybm, ps, plevs)

    # levels are pressure levels
    elif mv_plv.long_name.lower().find('pressure') != -1 or mv_plv.long_name.lower().find('isobaric') != -1:
        mv_p = utils.pressure_to_plevs(mv, plevs)

    else:
        raise RuntimeError(
            "Vertical level is neither hybrid nor pressure. Aborting.")
    
    return mv_p

def _regrid_and_plot(var, region, test_var, ref_var, ref_name, land_frac, ocean_frac, parameter):
    """
    Tentative implementation.
    For each region, regrid the test and ref data.
    Then take the difference and plot it.
    """
    print("Selected region: {}".format(region))

    test_var_domain, ref_var_domain = utils.select_region(
        region, test_var, ref_var, land_frac, ocean_frac, parameter)

    # Regrid towards lower resolution of two variables for
    # calculating difference
    test_var_reg, ref_var_reg = utils.regrid_to_lower_res(
        test_var_domain, ref_var_domain, parameter.regrid_tool, parameter.regrid_method)

    # A special case.
    if var == 'TREFHT_LAND' or var == 'SST':
        if ref_name == 'WILLMOTT':
            ref_var_reg = MV2.masked_where(
                ref_var_reg == ref_var_reg.fill_value, ref_var_reg)

        land_mask = MV2.logical_or(test_var_reg.mask, ref_var_reg.mask)
        test_var_reg = MV2.masked_where(land_mask, test_var_reg)
        ref_var_reg = MV2.masked_where(land_mask, ref_var_reg)

    # Plotting
    diff = test_var_reg - ref_var_reg
    metrics_dict = create_metrics(
        ref_var_domain, test_var_domain, ref_var_reg, test_var_reg, diff)

    metrics_dict['unit'] = test_var_reg.units

    fnm = os.path.join(get_output_dir(
        parameter.current_set, parameter), parameter.output_file)
    with open(fnm + '.json' , 'w') as outfile:
            json.dump(metrics_dict,outfile)
    print('Metrics saved in: ' + fnm + '.json')

    parameter.var_region = region
    plot(parameter.current_set, ref_var_domain,
            test_var_domain, diff, metrics_dict, parameter)
    utils.save_ncfiles(
        parameter.current_set, test_var_domain, ref_var_domain, diff, parameter)

def run_diag(parameter):
    parameter.reference_data_path
    parameter.test_data_path

    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = parameter.ref_name
    regions = parameter.regions

    for season in seasons:
        try:
            test_fnm = utils.get_test_filename(parameter, season)
            ref_fnm = utils.get_ref_filename(parameter, season)
        except IOError as e:
            print(e)
            traceback.print_exc()
            # the file for the current parameters wasn't found, move to next
            # parameters
            continue

        print('test file: {}'.format(test_fnm))
        print('reference file: {}'.format(ref_fnm))

        test_file = cdms2.open(test_fnm)
        ref_file = cdms2.open(ref_fnm)
        
        parameter.test_name_yrs = parameter.short_test_name if parameter.short_test_name else parameter.test_name

        try:
            yrs_averaged =  test_file.getglobal('yrs_averaged')
            parameter.test_name_yrs = parameter.test_name_yrs + ' (' + yrs_averaged +')'

        except:
            print('No yrs_averaged exists in global attributes')
            parameter.test_name_yrs = parameter.test_name_yrs

        land_frac = _get_land_or_ocean_frac(test_file, 'LANDFRAC')
        ocean_frac = _get_land_or_ocean_frac(test_file, 'OCNFRAC')

        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var
            test_var = acme.process_derived_var(
                var, acme.derived_variables, test_file, parameter)
            ref_var = acme.process_derived_var(
                var, acme.derived_variables, ref_file, parameter)

            parameter.viewer_descr[var] = test_var.long_name if hasattr(
                test_var, 'long_name') else 'No long_name attr in test data.'

            _handle_special_cases_for_ref(var, ref_name, season, ref_var)

            if test_var.getLevel() and ref_var.getLevel():  # for variables with z axis:
                plev = parameter.plevs

                print('Selected pressure level: {}'.format(plev))
                test_var_p = _convert_to_pressure_levels(test_var, test_file, plev)
                ref_var_p = _convert_to_pressure_levels(ref_var, ref_file, plev)

                # select plev
                for ilev in range(len(plev)):
                    test_var = test_var_p[ilev, ]
                    ref_var = ref_var_p[ilev, ]

                    for region in regions:
                        parameter.output_file = '-'.join(
                            [ref_name, var, str(int(plev[ilev])), season, region])
                        parameter.main_title = str(
                            ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))

                        _regrid_and_plot(var, region, test_var, ref_var, ref_name, land_frac, ocean_frac, parameter)

            # for variables without z axis:
            elif test_var.getLevel() is None and ref_var.getLevel() is None:

                for region in regions:
                    parameter.output_file = '-'.join(
                        [ref_name, var, str(int(plev[ilev])), season, region])
                    parameter.main_title = str(
                        ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))

                    _regrid_and_plot(var, region, test_var, ref_var, ref_name, land_frac, ocean_frac, parameter)

            else:
                raise RuntimeError(
                    "Dimensions of two variables are difference. Abort")

        ref_file.close()
        test_file.close()

    return parameter
