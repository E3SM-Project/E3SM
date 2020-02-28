from __future__ import print_function

import os
import numpy
import cdutil
import cdms2
import MV2
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean
from acme_diags.driver import utils
from acme_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter


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
        'mean': mean(ref, axis='yz')
    }
    metrics_dict['test'] = {
        'min': min_cdms(test),
        'max': max_cdms(test),
        'mean': mean(test, axis='yz')
    }

    metrics_dict['diff'] = {
        'min': min_cdms(diff),
        'max': max_cdms(diff),
        'mean': mean(diff, axis='yz')
    }
    metrics_dict['misc'] = {
        'rmse': rmse(test_regrid, ref_regrid, axis='yz'),
        'corr': corr(test_regrid, ref_regrid, axis='yz')
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
            land_frac = test_data.get_climo_variable('LANDFRAC', season)
            ocean_frac = test_data.get_climo_variable('OCNFRAC', season)
        except:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            with cdms2.open(mask_path) as f:
                land_frac = f('LANDFRAC')
                ocean_frac = f('OCNFRAC')

        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var

            mv1 = test_data.get_climo_variable(var, season)
            mv2 = ref_data.get_climo_variable(var, season)

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
                # Since the default is now stored in ZonalMean2dParameter,
                # we must get it from there if the plevs param is blank.
                plevs = parameter.plevs
                if (isinstance(plevs, numpy.ndarray) and not plevs.all()) or \
                    (not isinstance(plevs, numpy.ndarray) and not plevs):
                    plevs = ZonalMean2dParameter().plevs

                #print('Selected pressure level: {}'.format(plevs))

                mv1_p = utils.general.convert_to_pressure_levels(mv1, plevs, test_data, var, season)
                mv2_p = utils.general.convert_to_pressure_levels(mv2, plevs, ref_data, var, season)

                mv1_p = cdutil.averager(mv1_p, axis='x')
                mv2_p = cdutil.averager(mv2_p, axis='x')

                parameter.output_file = '-'.join(
                    [ref_name, var, season, parameter.regions[0]])
                parameter.main_title = str(
                    ' '.join([var, season]))
  

                # Regrid towards the lower resolution of the two
                # variables for calculating the difference.
                if len(mv1_p.getLatitude()) < len(mv2_p.getLatitude()):
                    mv1_reg = mv1_p
                    lev_out = mv1_p.getLevel()
                    lat_out = mv1_p.getLatitude()
                    mv2_reg = mv2_p.crossSectionRegrid(lev_out, lat_out)
                    # Apply the mask back, since crossSectionRegrid
                    # doesn't preserve the mask.
                    mv2_reg = MV2.masked_where(
                        mv2_reg == mv2_reg.fill_value, mv2_reg)
                elif len(mv1_p.getLatitude()) > len(mv2_p.getLatitude()):
                    mv2_reg = mv2_p
                    lev_out = mv2_p.getLevel()
                    lat_out = mv2_p.getLatitude()
                    mv1_reg = mv1_p.crossSectionRegrid(lev_out, lat_out)
                    # Apply the mask back, since crossSectionRegrid
                    # doesn't preserve the mask.
                    mv1_reg = MV2.masked_where(
                        mv1_reg == mv1_reg.fill_value, mv1_reg)
                else: 
                    mv1_reg = mv1_p
                    mv2_reg = mv2_p

                diff = mv1_reg - mv2_reg
                metrics_dict = create_metrics(
                    mv2_p, mv1_p, mv2_reg, mv1_reg, diff)

                parameter.var_region = 'global'

                plot(parameter.current_set, mv2_p, mv1_p,
                     diff, metrics_dict, parameter)
                utils.general.save_ncfiles(parameter.current_set,
                                   mv1_p, mv2_p, diff, parameter)

            # For variables without a z-axis.
            elif mv1.getLevel() is None and mv2.getLevel() is None:
                for region in regions:
                    #print("Selected region: {}".format(region))

                    mv1_domain = utils.general.select_region(region, mv1, land_frac, ocean_frac, parameter)
                    mv2_domain = utils.general.select_region(region, mv2, land_frac, ocean_frac, parameter)

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
                    parameter.var_region = region

                    plot(parameter.current_set, mv2_domain,
                         mv1_domain, diff, metrics_dict, parameter)
                    utils.general.save_ncfiles(parameter.current_set,
                                       mv1_domain, mv2_domain, diff, parameter)

            else:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting.")

    return parameter
