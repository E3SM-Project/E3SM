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
        
        if parameter.short_test_name:
            parameter.test_name_yrs = parameter.short_test_name
        else:
            parameter.test_name_yrs = parameter.test_name
        
        try:
            yrs_averaged =  test_file.getglobal('yrs_averaged')
            parameter.test_name_yrs = parameter.test_name_yrs + ' (' + yrs_averaged +')'

        except:
            print('No yrs_averaged exists in global attributes')
            parameter.test_name_yrs = parameter.test_name_yrs

        # save land/ocean fraction for masking
        try:
            land_frac = test_file('LANDFRAC')
            ocean_frac = test_file('OCNFRAC')
        except BaseException:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            f0 = cdms2.open(mask_path)
            land_frac = f0('LANDFRAC')
            ocean_frac = f0('OCNFRAC')
            f0.close()

        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var
            test_var = acme.process_derived_var(
                var, acme.derived_variables, test_file, parameter)
            ref_var = acme.process_derived_var(
                var, acme.derived_variables, ref_file, parameter)

            parameter.viewer_descr[var] = test_var.long_name if hasattr(
                test_var, 'long_name') else 'No long_name attr in test data.'

            # special case, cdms didn't properly convert mask with fill value
            # -999.0, filed issue with denise
            if ref_name == 'WARREN':
                # this is cdms2 for bad mask, denise's fix should fix
                ref_var = MV2.masked_where(ref_var == -0.9, ref_var)
                # following should move to derived variable
            if ref_name == 'AIRS':
                # ref_var = MV2.masked_where(ref_var==ref_var.fill_value, ref_var)
                # this is cdms2 for bad mask, denise's fix should fix
                ref_var = MV2.masked_where(ref_var > 1e+20, ref_var)
            if ref_name == 'WILLMOTT' or ref_name == 'CLOUDSAT':
                # ref_var = MV2.masked_where(ref_var==ref_var.fill_value, ref_var)
                # this is cdms2 for bad mask, denise's fix should fix
                ref_var = MV2.masked_where(ref_var == -999., ref_var)

                # following should move to derived variable
                if var == 'PRECT_LAND':
                    days_season = {'ANN': 365, 'DJF': 90,
                                   'MAM': 92, 'JJA': 92, 'SON': 91}
                    # test_var = test_var * days_season[season] * 0.1 #following AMWG
                    # approximate way to convert to seasonal cumulative
                    # precipitation, need to have solution in derived variable,
                    # unit convert from mm/day to cm
                    ref_var = ref_var / days_season[season] / \
                        0.1  # convert cm to mm/day instead
                    ref_var.units = 'mm/day'

            if test_var.getLevel() and ref_var.getLevel():  # for variables with z axis:
                plev = parameter.plevs
                print('Selected pressure level: {}'.format(plev))
                f_ins = [test_file, ref_file]
                for f_ind, mv in enumerate([test_var, ref_var]):
                    mv_plv = mv.getLevel()
                    # var(time,lev,lon,lat) convert from hybrid level to
                    # pressure
                    if mv_plv.long_name.lower().find('hybrid') != -1:
                        f_in = f_ins[f_ind]
                        hyam = f_in('hyam')
                        hybm = f_in('hybm')
                        ps = f_in('PS')  # Pa

                        mv_p = utils.hybrid_to_plevs(mv, hyam, hybm, ps, plev)

                    # levels are presure levels
                    elif mv_plv.long_name.lower().find('pressure') != -1 or mv_plv.long_name.lower().find('isobaric') != -1:
                        mv_p = utils.pressure_to_plevs(mv, plev)

                    else:
                        raise RuntimeError(
                            "Vertical level is neither hybrid nor pressure. Abort")
                    if f_ind == 0:
                        test_var_p = mv_p
                    if f_ind == 1:
                        ref_var_p = mv_p
                # select plev
                for ilev in range(len(plev)):
                    test_var = test_var_p[ilev, ]
                    ref_var = ref_var_p[ilev, ]

                    for region in regions:
                        print("Selected region: {}".format(region))

                        test_var_domain, ref_var_domain = utils.select_region(
                            region, test_var, ref_var, land_frac, ocean_frac, parameter)

                        parameter.output_file = '-'.join(
                            [ref_name, var, str(int(plev[ilev])), season, region])
                        parameter.main_title = str(
                            ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))

                        # Regrid towards lower resolution of two variables for
                        # calculating difference
                        test_var_reg, ref_var_reg = utils.regrid_to_lower_res(
                            test_var_domain, ref_var_domain, parameter.regrid_tool, parameter.regrid_method)

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

            # for variables without z axis:
            elif test_var.getLevel() is None and ref_var.getLevel() is None:

                # select region
                if len(regions) == 0:
                    regions = ['global']

                for region in regions:
                    print("Selected region: {}".format(region))

                    test_var_domain, ref_var_domain = utils.select_region(
                        region, test_var, ref_var, land_frac, ocean_frac, parameter)

                    parameter.output_file = '-'.join(
                        [ref_name, var, season, region])
                    parameter.main_title = str(' '.join([var, season, region]))

                    # regrid towards lower resolution of two variables for
                    # calculating difference
                    test_var_reg, ref_var_reg = utils.regrid_to_lower_res(
                        test_var_domain, ref_var_domain, parameter.regrid_tool, parameter.regrid_method)

                    # if var is 'SST' or var is 'TREFHT_LAND': #special case

                    if var == 'TREFHT_LAND'or var == 'SST':  # use "==" instead of "is"
                        if ref_name == 'WILLMOTT':
                            ref_var_reg = MV2.masked_where(
                                ref_var_reg == ref_var_reg.fill_value, ref_var_reg)
                            print(ref_name)

                            # if mv.mask is False:
                            #    mv = MV2.masked_less_equal(mv, mv._FillValue)
                            #    print("*************",mv.count())
                        land_mask = MV2.logical_or(test_var_reg.mask, ref_var_reg.mask)
                        test_var_reg = MV2.masked_where(land_mask, test_var_reg)
                        ref_var_reg = MV2.masked_where(land_mask, ref_var_reg)

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
                    utils.save_ncfiles(parameter.current_set,
                                       test_var_domain, ref_var_domain, diff, parameter)

            else:
                raise RuntimeError(
                    "Dimensions of two variables are difference. Abort")
        ref_file.close()
        test_file.close()
    return parameter
    