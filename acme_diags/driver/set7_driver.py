#!/usr/bin/env python
import json
import copy
import glob
import os
import sys
import fnmatch
import numpy
import cdutil
import genutil
import cdms2
import MV2
from acme_diags.acme_viewer import create_viewer
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.derivations.default_regions import regions_specs
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean

def findfile(path_name, data_name, season):
    """locate file name based on data_name and season"""
    dir_files = os.listdir(path_name)
    for filename in dir_files:
        if filename.startswith(data_name) and season in filename:
            return path_name+filename
    raise IOError(
        "No file found based on given path_name and data_name")

def regrid_to_lower_res(mv1, mv2, regrid_tool, regrid_method):
    """regrid transient variable toward lower resolution of two variables"""

    axes1 = mv1.getAxisList()
    axes2 = mv2.getAxisList()

    # use nlat to decide data resolution, higher number means higher data
    # resolution. For the difference plot, regrid toward lower resolution
    if len(axes1[1]) <= len(axes2[1]):
        mv_grid = mv1.getGrid()
        mv1_reg = mv1
        mv2_reg = mv2.regrid(mv_grid, regridTool=regrid_tool,
                             regridMethod=regrid_method)
    else:
        mv_grid = mv2.getGrid()
        mv2_reg = mv2
        mv1_reg = mv1.regrid(mv_grid, regridTool=regrid_tool,
                             regridMethod=regrid_method)
    return mv1_reg, mv2_reg


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """ Creates the mean, max, min, rmse, corr in a dictionary """
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


def mask_by(input_var, maskvar, low_limit=None, high_limit=None):
    """masks a variable var to be missing except where maskvar>=low_limit and maskvar<=high_limit. 
    None means to omit the constrint, i.e. low_limit = -infinity or high_limit = infinity. var is changed and returned; we don't make a new variable.
    var and maskvar: dimensioned the same variables.
    low_limit and high_limit: scalars.
    """
    var = copy.deepcopy(input_var)
    if low_limit is None and high_limit is None:
        return var
    if low_limit is None and high_limit is not None:
        maskvarmask = maskvar > high_limit
    elif low_limit is not None and high_limit is None:
        maskvarmask = maskvar < low_limit
    else:
        maskvarmask = (maskvar < low_limit) | (maskvar > high_limit)
    if var.mask is False:
        newmask = maskvarmask
    else:
        newmask = var.mask | maskvarmask
    var.mask = newmask
    return var


def save_ncfiles(test, ref, diff, parameter):
    ''' Saves the test, reference, and difference nc files. '''
    # Save files being plotted
    # Set cdms preferences - no compression, no shuffling, no complaining
    cdms2.setNetcdfDeflateFlag(1)
    # 1-9, min to max - Comes at heavy IO (read/write time cost)
    cdms2.setNetcdfDeflateLevelFlag(0)
    cdms2.setNetcdfShuffleFlag(0)
    cdms2.setCompressionWarnings(0)  # Turn off warning messages
    # Save test file
    file_test = cdms2.open(parameter.results_dir + '/' + parameter.case_id + 
                           '/' + parameter.output_file + '_test.nc', 'w+')
    test.id = parameter.var_id
    file_test.write(test)
    file_test.close()

    # Save reference file
    file_ref = cdms2.open(parameter.results_dir + '/' + parameter.case_id + 
                           '/' + parameter.output_file + '_ref.nc', 'w+')
    ref.id = parameter.var_id
    file_ref.write(ref)
    file_ref.close()

    # Save difference file
    file_diff = cdms2.open(parameter.results_dir + '/' + parameter.case_id + 
                           '/' + parameter.output_file + '_diff.nc', 'w+')
    diff.id = parameter.var_id + '(test - reference)'
    file_diff.write(diff)
    file_diff.close()

def run_diag(parameter):
    if not os.path.exists(os.path.join(parameter.results_dir, parameter.case_id)):
        os.makedirs(os.path.join(parameter.results_dir, parameter.case_id))

    reference_data_path = parameter.reference_data_path
    test_data_path = parameter.test_data_path

    variables = parameter.variables
    seasons = parameter.season
    ref_name = parameter.ref_name
    test_name = parameter.test_name
    regions = parameter.region

    for season in seasons:
        if hasattr(parameter, 'test_path'):
            filename1 = parameter.test_path
            print filename1
        else:
            try:
                filename1 = findfile(test_data_path, test_name, season)
                print filename1
            except IOError:
                print('No file found for {} and {}'.format(test_name, season))
                continue

        if hasattr(parameter, 'reference_path'):
            filename2 = parameter.reference_path
            print filename2
        else:
            try:
                filename2 = findfile(reference_data_path, ref_name, season)
                print filename2
            except IOError:
                print('No file found for {} and {}'.format(test_name, season))
                continue

        f_mod = cdms2.open(filename1)
        f_obs = cdms2.open(filename2)

        for var in variables: 
            print '***********', variables
            print var
            parameter.var_id = var
            mv1 = acme.process_derived_var(
                var, acme.derived_variables, f_mod, parameter)
            mv2 = acme.process_derived_var(
                var, acme.derived_variables, f_obs, parameter)
    
            # special case, cdms didn't properly convert mask with fill value
            # -999.0, filed issue with denise
            if ref_name == 'WARREN':
                # this is cdms2 for bad mask, denise's fix should fix
                mv2 = MV2.masked_where(mv2 == -0.9, mv2)
                # following should move to derived variable
            if ref_name == 'WILLMOTT' or ref_name == 'CLOUDSAT':
                print mv2.fill_value
                # mv2=MV2.masked_where(mv2==mv2.fill_value,mv2)
                # this is cdms2 for bad mask, denise's fix should fix
                mv2 = MV2.masked_where(mv2 == -999., mv2)
    
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
    
            # for variables without z axis:
            if mv1.getLevel() == None and mv2.getLevel() == None:
                if len(regions) == 0:
                    regions = ['global']
    
                for region in regions:
                    print region
                    domain = None
                    # if region != 'global':
                    if region.find('land') != -1 or region.find('ocean') != -1:
                        if region.find('land') != -1:
                            land_ocean_frac = f_mod('LANDFRAC')
                        elif region.find('ocean') != -1:
                            land_ocean_frac = f_mod('OCNFRAC')
                        region_value = regions_specs[region]['value']
                        print 'region_value', region_value, mv1
    
                        mv1_domain = mask_by(
                            mv1, land_ocean_frac, low_limit=region_value)
                        mv2_domain = mv2.regrid(
                            mv1.getGrid(), parameter.regrid_tool, parameter.regrid_method)
                        mv2_domain = mask_by(
                            mv2_domain, land_ocean_frac, low_limit=region_value)
                    else:
                        mv1_domain = mv1
                        mv2_domain = mv2
    
                    print region
                    try:
                        domain = regions_specs[region]['domain']
                        print domain
                    except:
                        print("no domain selector")
                    mv1_domain = mv1_domain(domain)
                    mv2_domain = mv2_domain(domain)
                    mv1_domain.units = mv1.units
                    mv2_domain.units = mv1.units
    
                    parameter.output_file = '-'.join(
                        [ref_name, var, season, region])
                    parameter.main_title = str(' '.join([var, season, region]))
    
                    # regrid towards lower resolution of two variables for
                    # calculating difference
                    mv1_reg, mv2_reg = regrid_to_lower_res(
                        mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)
    
                    # if var is 'SST' or var is 'TREFHT_LAND': #special case
    
                    if var == 'TREFHT_LAND'or var == 'SST':  # use "==" instead of "is"
                        if ref_name == 'WILLMOTT':
                            mv2_reg = MV2.masked_where(
                                mv2_reg == mv2_reg.fill_value, mv2_reg)
                            print ref_name
    
                            # if mv.mask is False:
                            #    mv = MV2.masked_less_equal(mv, mv._FillValue)
                            #    print "*************",mv.count()
                        land_mask = MV2.logical_or(mv1_reg.mask, mv2_reg.mask)
                        mv1_reg = MV2.masked_where(land_mask, mv1_reg)
                        mv2_reg = MV2.masked_where(land_mask, mv2_reg)
    
                    diff = mv1_reg - mv2_reg
                    metrics_dict = create_metrics(
                        mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)
                    parameter.var_region = region
                    plot('7', mv2_domain, mv1_domain, diff, metrics_dict, parameter)

                    save_ncfiles(mv1_domain, mv2_domain, diff, parameter)
    
            elif mv1.getLevel() and mv2.getLevel():  # for variables with z axis:
                plev = parameter.plevs
                print 'selected pressure level', plev
                f_mod = cdms2.open(filename1)
                for filename in [filename1, filename2]:
                    f_in = cdms2.open(filename)
                    mv = f_in[var]  # Square brackets for metadata preview
                    mv_plv = mv.getLevel()
    
                    # var(time,lev,lon,lat) convert from hybrid level to pressure
                    if mv_plv.long_name.lower().find('hybrid') != -1:
                        # Parentheses actually load the transient variable
                        mv = f_in(var)
                        hyam = f_in('hyam')
                        hybm = f_in('hybm')
                        ps = f_in('PS') / 100.  # convert unit from 'Pa' to mb
                        p0 = 1000.  # mb
                        levels_orig = cdutil.vertical.reconstructPressureFromHybrid(
                            ps, hyam, hybm, p0)
                        levels_orig.units = 'mb'
                        mv_p = cdutil.vertical.logLinearInterpolation(
                            mv, levels_orig, plev)
    
                    elif mv_plv.long_name.lower().find('pressure') != -1 or mv_plv.long_name.lower().find('isobaric') != -1:  # levels are presure levels
    
                        mv = f_in(var)
                        #mv = acme.process_derived_var(var, acme.derived_variables, f_in)
                        if ref_name == 'AIRS':
                            # mv2=MV2.masked_where(mv2==mv2.fill_value,mv2)
                            print mv.fill_value
                            # this is cdms2 for bad mask, denise's fix should fix
                            mv = MV2.masked_where(mv == mv.fill_value, mv)
                        # Construct pressure level for interpolation
                        levels_orig = MV2.array(mv_plv[:])
                        levels_orig.setAxis(0, mv_plv)
                        # grow 1d levels_orig to mv dimention
                        mv, levels_orig = genutil.grower(mv, levels_orig)
                        # levels_orig.info()
                        # logLinearInterpolation only takes positive down plevel:
                        # "I :      interpolation field (usually Pressure or depth)
                        # from TOP (level 0) to BOTTOM (last level), i.e P value
                        # going up with each level"
                        mv_p = cdutil.vertical.logLinearInterpolation(
                            mv[:, ::-1, ], levels_orig[:, ::-1, ], plev)
    
                    else:
                        raise RuntimeError(
                            "Vertical level is neither hybrid nor pressure. Abort")
    
                    if filename == filename1:
                        mv1_p = mv_p
    
                    if filename == filename2:
                        mv2_p = mv_p
    
                for ilev in range(len(plev)):
                    mv1 = mv1_p[:, ilev, ]
                    mv2 = mv2_p[:, ilev, ]
                    if len(regions) == 0:
                        regions = ['global']
    
                    for region in regions:
                        print region
                        domain = None
                        # if region != 'global':
                        if region.find('land') != -1 or region.find('ocean') != -1:
                            if region.find('land') != -1:
                                land_ocean_frac = f_mod('LANDFRAC')
                            elif region.find('ocean') != -1:
                                land_ocean_frac = f_mod('OCNFRAC')
                            region_value = regions_specs[region]['value']
                            print 'region_value', region_value, mv1
    
                            mv1_domain = mask_by(
                                mv1, land_ocean_frac, low_limit=region_value)
                            mv2_domain = mv2.regrid(
                                mv1.getGrid(), parameter.regrid_tool, parameter.regrid_method)
                            mv2_domain = mask_by(
                                mv2_domain, land_ocean_frac, low_limit=region_value)
                        else:
                            mv1_domain = mv1
                            mv2_domain = mv2
    
                        print region
                        try:
                            # if region.find('global') == -1:
                            domain = regions_specs[region]['domain']
                            print domain
                        except:
                            print ("no domain selector")
                        mv1_domain = mv1_domain(domain)
                        mv2_domain = mv2_domain(domain)
                        mv1_domain.units = mv1.units
                        mv2_domain.units = mv1.units
    
                        parameter.output_file = '-'.join(
                            [ref_name, var, str(int(plev[ilev])), season, region])
                        parameter.main_title = str(
                            ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))
    
                        # Regrid towards lower resolution of two variables for
                        # calculating difference
                        mv1_reg, mv2_reg = regrid_to_lower_res(
                            mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)
    
                        # Plotting
                        diff = mv1_reg - mv2_reg
                        metrics_dict = create_metrics(
                            mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)

                        parameter.var_region = region
                        plot('7', mv2_domain, mv1_domain, diff, metrics_dict, parameter)
                        save_ncfiles(mv1_domain, mv2_domain, diff, parameter)

                f_in.close()
            
            else:
                raise RuntimeError(
                    "Dimensions of two variables are difference. Abort")
        f_obs.close()
        f_mod.close()
