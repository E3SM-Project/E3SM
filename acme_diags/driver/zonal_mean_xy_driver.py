from __future__ import print_function

import os
import sys
import cdms2
import MV2
import cdutil
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean
from acme_diags.driver import utils, dataset


def regrid_to_lower_res_1d(mv1, mv2):
    """Regrid 1-D transient variable toward lower resolution of two variables."""
    import numpy

    if mv1 is None or mv2 is None:
        return None
    missing = mv1.get_fill_value()
    axis1 = mv1.getAxisList()[0]
    axis2 = mv2.getAxisList()[0]
    if len(axis1) <= len(axis2):
        mv1_reg = mv1
        b0 = numpy.interp(axis1[:], axis2[:], mv2[:],
                          left=missing, right=missing)
        b0_mask = numpy.interp(
            axis1[:], axis2[:], mv2.mask[:], left=missing, right=missing)
        mv2_reg = cdms2.createVariable(b0, mask=[True if bb == missing or bb_mask != 0.0 else False for (
            bb, bb_mask) in zip(b0[:], b0_mask[:])], axes=[axis1])
    else:
        a0 = numpy.interp(axis2[:], axis1[:], mv1[:],
                          left=missing, right=missing)
        a0_mask = numpy.interp(
            axis2[:], axis1[:], mv1.mask[:], left=missing, right=missing)
        mv1_reg = cdms2.createVariable(a0, mask=[True if aa == missing or aa_mask != 0.0 else False for (
            aa, aa_mask) in zip(a0[:], a0_mask[:])], axes=[axis2])
        mv2_reg = mv2

    return mv1_reg, mv2_reg


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


def _convert_to_pressure_levels(mv, plevs, dataset, season):
    """
    Given either test or reference data with a z-axis,
    convert to the desired pressure levels.
    """
    mv_plv = mv.getLevel()
    # var(time,lev,lon,lat) convert from hybrid level to pressure
    if mv_plv.long_name.lower().find('hybrid') != -1:
        hyam = dataset.get_variable('hyam', season)
        hybm = dataset.get_variable('hybm', season)
        ps = dataset.get_variable('PS', season)  # Pa
        mv_p = utils.hybrid_to_plevs(mv, hyam, hybm, ps, plevs)

    # levels are pressure levels
    elif mv_plv.long_name.lower().find('pressure') != -1 or mv_plv.long_name.lower().find('isobaric') != -1:
        mv_p = utils.pressure_to_plevs(mv, plevs)

    else:
        raise RuntimeError(
            "Vertical level is neither hybrid nor pressure. Aborting.")

    return mv_p


def run_diag(parameter):
    parameter.reference_data_path
    parameter.test_data_path

    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = parameter.ref_name
    regions = parameter.regions

    test_data = dataset.Dataset(parameter, test=True)
    ref_data = dataset.Dataset(parameter, ref=True)    
    # test_data.get_variable(var, season)

    for season in seasons:
        if parameter.short_test_name:
            parameter.test_name_yrs = parameter.short_test_name
        else:
            parameter.test_name_yrs = parameter.test_name

        
        try:
            # yrs_averaged = f_mod.getglobal('yrs_averaged')
            # parameter.test_name_yrs = parameter.test_name_yrs + ' (' + yrs_averaged +')'
            # Raising an exception b/c we don't have an f_mod variable.
            raise Exception()
        except:
            print('No yrs_averaged exists in global attributes')
            parameter.test_name_yrs = parameter.test_name_yrs

        # save land/ocean fraction for masking
        try:
            test_data.get_variable('LANDFRAC', season)
            test_data.get_variable('OCNFRAC', season)
        except:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            with cdms2.open(mask_path) as f:
                f('LANDFRAC')
                f('OCNFRAC')

        
        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var
            mv1 = test_data.get_variable(var, season)
            mv2 = ref_data.get_variable(var, season)

            parameter.viewer_descr[var] = mv1.long_name if hasattr(
                mv1, 'long_name') else 'No long_name attr in test data.'

            # special case, cdms didn't properly convert mask with fill value
            # -999.0, filed issue with denise
            if ref_name == 'WARREN':
                # this is cdms2 for bad mask, denise's fix should fix
                mv2 = MV2.masked_where(mv2 == -0.9, mv2)
                # following should move to derived variable
            if ref_name == 'AIRS':
                # mv2=MV2.masked_where(mv2==mv2.fill_value,mv2)
                # this is cdms2 for bad mask, denise's fix should fix
                mv2 = MV2.masked_where(mv2 > 1e+20, mv2)
            if ref_name == 'WILLMOTT' or ref_name == 'CLOUDSAT':
                print(mv2.fill_value)
                # mv2=MV2.masked_where(mv2==mv2.fill_value,mv2)
                # this is cdms2 for bad mask, denise's fix should fix
                mv2 = MV2.masked_where(mv2 == -999., mv2)
                print(mv2.fill_value)

                # following should move to derived variable
                if var == 'PRECT_LAND':
                    days_season = {'ANN': 365, 'DJF': 90,
                                   'MAM': 92, 'JJA': 92, 'SON': 91}
                    # mv1 = mv1 * days_season[season] * 0.1 #following AMWG
                    # approximate way to convert to seasonal cumulative
                    # precipitation, need to have solution in derived variable,
                    # unit convert from mm/day to cm
                    mv2 = mv2 / days_season[season] / \
                        0.1  # convert cm to mm/day instead
                    mv2.units = 'mm/day'

            if mv1.getLevel() and mv2.getLevel():  # for variables with z axis:
                plev = parameter.plevs
                print('Selected pressure level: {}'.format(plev))
                
                mv1_p = _convert_to_pressure_levels(mv1, plev, test_data, season)
                mv2_p = _convert_to_pressure_levels(mv2, plev, test_data, season)

                # select plev
                for ilev in range(len(plev)):
                    mv1 = mv1_p[ilev, ]
                    mv2 = mv2_p[ilev, ]

                    # select region
                    if len(regions) == 0:
                        regions = ['global']

                    for region in regions:
                        print("Selected region: {}".format(region))
                        mv1_zonal = cdutil.averager(mv1, axis='x')
                        mv2_zonal = cdutil.averager(mv2, axis='x')

                        # Regrid towards lower resolution of two variables for
                        # calculating difference
                        print(mv1_zonal.shape, mv2_zonal.shape)
                        mv1_reg, mv2_reg = regrid_to_lower_res_1d(
                            mv1_zonal, mv2_zonal)

                        diff = mv1_reg - mv2_reg
                        parameter.output_file = '-'.join(
                            [ref_name, var, str(int(plev[ilev])), season, region])
                        parameter.main_title = str(
                            ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))

                        parameter.var_region = region
                        plot(parameter.current_set, mv2_zonal,
                             mv1_zonal, diff, {}, parameter)
                        utils.save_ncfiles(
                            parameter.current_set, mv1_zonal, mv2_zonal, diff, parameter)

#                        mv1_domain, mv2_domain = utils.select_region(region, mv1, mv2, land_frac,ocean_frac,parameter)
#
#                        parameter.output_file = '-'.join(
#                            [ref_name, var, str(int(plev[ilev])), season, region])
#                        parameter.main_title = str(
#                            ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))
#
#                        # Regrid towards lower resolution of two variables for
#                        # calculating difference
#                        mv1_reg, mv2_reg = utils.regrid_to_lower_res(
#                            mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)
#
#                        # Plotting
#                        diff = mv1_reg - mv2_reg
#                        metrics_dict = create_metrics(
#                            mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)
#
#                        parameter.var_region = region
#                        plot(parameter.current_set, mv2_zonal, mv1_zonal, diff parameter)
#                        utils.save_ncfiles(parameter.current_set, mv1_domain, mv2_domain, diff, parameter)

            # for variables without z axis:
            elif mv1.getLevel() is None and mv2.getLevel() is None:

                # select region
                if len(regions) == 0:
                    regions = ['global']

                for region in regions:
                    print("Selected region: {}".format(region))
                    mv1_zonal = cdutil.averager(mv1, axis='x')
                    mv2_zonal = cdutil.averager(mv2, axis='x')

                    print(mv1_zonal.shape, mv2_zonal.shape)
                    mv1_reg, mv2_reg = regrid_to_lower_res_1d(
                        mv1_zonal, mv2_zonal)

                    diff = mv1_reg - mv2_reg

                    parameter.output_file = '-'.join(
                        [ref_name, var, season, region])
                    parameter.main_title = str(' '.join([var, season, region]))

                    parameter.var_region = region

                    plot(parameter.current_set, mv2_zonal,
                         mv1_zonal, diff, {}, parameter)
                    utils.save_ncfiles(parameter.current_set,
                                       mv1_zonal, mv2_zonal, diff, parameter)


#                    mv1_domain, mv2_domain = utils.select_region(region, mv1, mv2, land_frac,ocean_frac,parameter)
#
#                    parameter.output_file = '-'.join(
#                        [ref_name, var, season, region])
#                    parameter.main_title = str(' '.join([var, season, region]))
#
#                    # regrid towards lower resolution of two variables for
#                    # calculating difference
#                    mv1_reg, mv2_reg = utils.regrid_to_lower_res(
#                        mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)
#
#                    # if var is 'SST' or var is 'TREFHT_LAND': #special case
#
#                    if var == 'TREFHT_LAND'or var == 'SST':  # use "==" instead of "is"
#                        if ref_name == 'WILLMOTT':
#                            mv2_reg = MV2.masked_where(
#                                mv2_reg == mv2_reg.fill_value, mv2_reg)
#                            print(ref_name)
#
#                            # if mv.mask is False:
#                            #    mv = MV2.masked_less_equal(mv, mv._FillValue)
#                            #    print("*************",mv.count())
#                        land_mask = MV2.logical_or(mv1_reg.mask, mv2_reg.mask)
#                        mv1_reg = MV2.masked_where(land_mask, mv1_reg)
#                        mv2_reg = MV2.masked_where(land_mask, mv2_reg)
#
#                    diff = mv1_reg - mv2_reg
#                    metrics_dict = create_metrics(
#                        mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)
#                    parameter.var_region = region
#                    plot('5', mv2_domain, mv1_domain, diff, metrics_dict, parameter)
#                    utils.save_ncfiles('5', mv1_domain, mv2_domain, diff, parameter)

            else:
                raise RuntimeError(
                    "Dimensions of two variables are difference. Abort")

    return parameter
