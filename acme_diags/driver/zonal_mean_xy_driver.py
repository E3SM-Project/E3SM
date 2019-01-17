from __future__ import print_function

import os
import numpy
import cdms2
import MV2
import cdutil
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean
from acme_diags.driver import utils


def regrid_to_lower_res_1d(mv1, mv2):
    """Regrid 1-D transient variable toward lower resolution of two variables."""

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
                mv2_p = utils.general.convert_to_pressure_levels(mv2, plev, test_data, var, season)

                # Select plev.
                for ilev in range(len(plev)):
                    mv1 = mv1_p[ilev, ]
                    mv2 = mv2_p[ilev, ]

                    for region in regions:
                        print("Selected region: {}".format(region))
                        mv1_zonal = cdutil.averager(mv1, axis='x')
                        mv2_zonal = cdutil.averager(mv2, axis='x')

                        # Regrid towards the lower resolution of the two
                        # variables for calculating the difference.
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
                        utils.general.save_ncfiles(
                            parameter.current_set, mv1_zonal, mv2_zonal, diff, parameter)

            # For variables without a z-axis.
            elif mv1.getLevel() is None and mv2.getLevel() is None:
                for region in regions:
                    print("Selected region: {}".format(region))
                    mv1_zonal = cdutil.averager(mv1, axis='x')
                    mv2_zonal = cdutil.averager(mv2, axis='x')

                    mv1_reg, mv2_reg = regrid_to_lower_res_1d(
                        mv1_zonal, mv2_zonal)

                    diff = mv1_reg - mv2_reg

                    parameter.output_file = '-'.join(
                        [ref_name, var, season, region])
                    parameter.main_title = str(' '.join([var, season, region]))

                    parameter.var_region = region

                    plot(parameter.current_set, mv2_zonal,
                         mv1_zonal, diff, {}, parameter)
                    utils.general.save_ncfiles(parameter.current_set,
                                       mv1_zonal, mv2_zonal, diff, parameter)

            else:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting.")

    return parameter
