from __future__ import print_function

import os
import sys
import numpy
import cdutil
import cdms2
import MV2
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean
from acme_diags.driver import utils


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    orig_bounds = cdms2.getAutoBounds()
    cdms2.setAutoBounds(1)
    lev = ref.getLevel()
    if lev is not None:
        lev.setBounds(None)

    lev = test.getLevel()
    if lev is not None:
        lev.setBounds(None)

    lev = test_regrid.getLevel()
    if lev is not None:
        lev.setBounds(None)

    lev = ref_regrid.getLevel()
    if lev is not None:
        lev.setBounds(None)

    lev = diff.getLevel()
    if lev is not None:
        lev.setBounds(None)
    cdms2.setAutoBounds(orig_bounds)

    metrics_dict = {}
    metrics_dict['ref'] = {
        'min': min_cdms(ref),
        'max': max_cdms(ref),
        # 'mean': numpy.nan #mean(ref, axis='yz')
        'mean': mean(ref, axis='yz')
    }
    metrics_dict['test'] = {
        'min': min_cdms(test),
        'max': max_cdms(test),
        # 'mean': numpy.nan #mean(test, axis='yz')
        'mean': mean(test, axis='yz')
    }

    metrics_dict['diff'] = {
        'min': min_cdms(diff),
        'max': max_cdms(diff),
        # 'mean': numpy.nan #mean(diff, axis='yz')
        'mean': mean(diff, axis='yz')
    }
    metrics_dict['misc'] = {
        'rmse': rmse(test_regrid, ref_regrid, axis='yz'),
        'corr': corr(test_regrid, ref_regrid, axis='yz')
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
            filename1 = utils.get_test_filename(parameter, season)
            filename2 = utils.get_ref_filename(parameter, season)
        except IOError as e:
            print(e)
            # the file for the current parameters wasn't found, move to next
            # parameters
            continue

        print('test file: {}'.format(filename1))
        print('reference file: {}'.format(filename2))

        f_mod = cdms2.open(filename1)
        f_obs = cdms2.open(filename2)

        if parameter.short_test_name:
            parameter.test_name_yrs = parameter.short_test_name
        else:
            parameter.test_name_yrs = parameter.test_name

        try:
            yrs_averaged =  f_mod.getglobal('yrs_averaged')
            parameter.test_name_yrs = parameter.test_name_yrs + ' (' + yrs_averaged +')'

        except:
            print('No yrs_averaged exists in global attributes')
            parameter.test_name_yrs = parameter.test_name_yrs

        # save land/ocean fraction for masking
        try:
            land_frac = f_mod('LANDFRAC')
            ocean_frac = f_mod('OCNFRAC')
        except BaseException:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            f0 = cdms2.open(mask_path)
            land_frac = f0('LANDFRAC')
            ocean_frac = f0('OCNFRAC')
            f0.close()

        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var
            mv1 = acme.process_derived_var(
                var, acme.derived_variables, f_mod, parameter)
            mv2 = acme.process_derived_var(
                var, acme.derived_variables, f_obs, parameter)

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
                print(mv2.fill_value)
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
                if parameter.zonal_mean_2d_plevs:
                    plev = parameter.zonal_mean_2d_plevs
                else:
                    plev = numpy.logspace(2.0, 3.0, num=17)
                # plev = [30.,50.,70.,100.,150.,200.,250.,300.,400.,500.,600.,700.,775.,850.,925.,1000.]
                print('Selected pressure level: {}'.format(plev))
                f_ins = [f_mod, f_obs]
                for f_ind, mv in enumerate([mv1, mv2]):
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

                    # calculate zonal mean
                    mv_p = cdutil.averager(mv_p, axis='x')
                    if f_ind == 0:
                        mv1_p = mv_p
                    elif f_ind == 1:
                        mv2_p = mv_p
                parameter.output_file = '-'.join(
                    [ref_name, var, season, parameter.regions[0]])
                parameter.main_title = str(
                    ' '.join([var, season]))

                # Regrid towards lower resolution of two variables for
                # calculating difference
                if len(mv1_p.getLatitude()) <= len(mv2_p.getLatitude()):
                    mv1_reg = mv1_p
                    lev_out = mv1_p.getLevel()
                    lat_out = mv1_p.getLatitude()
                    mv2_reg = mv2_p.crossSectionRegrid(lev_out, lat_out)
                    # apply mask back, since crossSectionRegrid doesn't
                    # preserve mask
                    mv2_reg = MV2.masked_where(
                        mv2_reg == mv2_reg.fill_value, mv2_reg)
                    print(mv2_reg.fill_value)
                else:
                    mv2_reg = mv2_p
                    lev_out = mv2_p.getLevel()
                    lat_out = mv2_p.getLatitude()
                    mv1_reg = mv1_p.crossSectionRegrid(lev_out, lat_out)
                    # apply mask back, since crossSectionRegrid doesn't
                    # preserve mask
                    mv1_reg = MV2.masked_where(
                        mv1_reg == mv1_reg.fill_value, mv1_reg)

                # print(mv1_p.shape, mv2_p.shape)
                # mv1_reg, mv2_reg = utils.sregrid_to_lower_res(
                #    mv1_p, mv2_p, parameter.regrid_tool, parameter.regrid_method)
                # reg_mask = MV2.logical_and(mv1_reg.mask, mv2_reg.mask)
                # print('reg_mask', reg_mask[:,-1])
                # mv1_reg = MV2.masked_where(reg_mask,mv1_reg)
                # mv2_reg = MV2.masked_where(reg_mask,mv2_reg)
                # print(mv1_reg.shape)
                # print(mv2_reg.shape)
                # print(mv1_p[:,-1].mask)
                # print(mv2_p[:,-1].mask)

                # print(mv1_reg[:,-1].mask)
                # print(mv2_reg[:,-1].mask)
                # Plotting
                diff = mv1_reg - mv2_reg
                metrics_dict = create_metrics(
                    mv2_p, mv1_p, mv2_reg, mv1_reg, diff)

                parameter.var_region = 'global'

                plot(parameter.current_set, mv2_p, mv1_p,
                     diff, metrics_dict, parameter)
                utils.save_ncfiles(parameter.current_set,
                                   mv1_p, mv2_p, diff, parameter)

            # for variables without z axis:
            elif mv1.getLevel() is None and mv2.getLevel() is None:

                # select region
                if len(regions) == 0:
                    regions = ['global']

                for region in regions:
                    print("Selected region: {}".format(region))

                    mv1_domain, mv2_domain = utils.select_region(
                        region, mv1, mv2, land_frac, ocean_frac, parameter)

                    parameter.output_file = '-'.join(
                        [ref_name, var, season, region])
                    parameter.main_title = str(' '.join([var, season, region]))

                    # regrid towards lower resolution of two variables for
                    # calculating difference
                    mv1_reg, mv2_reg = utils.regrid_to_lower_res(
                        mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                    # if var is 'SST' or var is 'TREFHT_LAND': #special case

                    if var == 'TREFHT_LAND'or var == 'SST':  # use "==" instead of "is"
                        if ref_name == 'WILLMOTT':
                            mv2_reg = MV2.masked_where(
                                mv2_reg == mv2_reg.fill_value, mv2_reg)
                            print(ref_name)

                            # if mv.mask is False:
                            #    mv = MV2.masked_less_equal(mv, mv._FillValue)
                            #    print("*************",mv.count())
                        land_mask = MV2.logical_or(mv1_reg.mask, mv2_reg.mask)
                        mv1_reg = MV2.masked_where(land_mask, mv1_reg)
                        mv2_reg = MV2.masked_where(land_mask, mv2_reg)

                    diff = mv1_reg - mv2_reg
                    metrics_dict = create_metrics(
                        mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)
                    parameter.var_region = region

                    plot(parameter.current_set, mv2_domain,
                         mv1_domain, diff, metrics_dict, parameter)
                    utils.save_ncfiles(parameter.current_set,
                                       mv1_domain, mv2_domain, diff, parameter)

            else:
                raise RuntimeError(
                    "Dimensions of two variables are difference. Abort")
        f_obs.close()
        f_mod.close()
    return parameter
