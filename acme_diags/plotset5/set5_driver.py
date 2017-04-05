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
from cdp.cdp_viewer import OutputViewer
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_parameter import ACMEParameter
from acme_diags.plotting.set5.plot import plot
from acme_diags.derivations import acme
from acme_diags.derivations.default_regions import regions_specs
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean


def make_parameters(orginal_parameter):
    """ Create multiple parameters given a list of
    parameters in a json and an original parameter """

    pth = os.path.join(sys.prefix, 'share', 'acme_diags', 'set5', 'set5_diags_AMWG_default.json')
    f_data = open(pth).read()
    json_file = json.loads(f_data)

    # add the custom_diags to the main default diags file
    if hasattr(original_parameter, 'custom_diags'):
        f_data = open(original_parameter.custom_diags).read()
        custom_json_data = json.loads(f_data)

        for key in custom_json_data:
            for single_run in custom_json_data[key]:
                json_file['set5'].append(single_run)

    parameters = []
    for key in json_file:
        for single_run in json_file[key]:
            p = ACMEParameter()
            for attr_name in single_run:
                setattr(p, attr_name, single_run[attr_name])

            # Add attributes of original_parameter to p
            p.__dict__.update(orginal_parameter.__dict__)
            parameters.append(p)
    return parameters


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


def add_page_and_top_row(viewer, parameters):
    """ Setup for OutputViewer """
    col_labels = ['Description']
    seasons = []
    for p in parameters:
        for s in p.season:
            if s not in seasons:
                seasons.append(s)
    for s in seasons:
        col_labels.append(s)

    viewer.add_page("Set 5", col_labels)


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
    file_test = cdms2.open(parameter.case_id + '/' +
                           parameter.output_file + '_test.nc', 'w+')
    test.id = parameter.variables
    file_test.write(test)
    file_test.close()

    # Save reference file
    file_ref = cdms2.open(parameter.case_id + '/' +
                          parameter.output_file + '_ref.nc', 'w+')
    ref.id = parameter.variables
    file_ref.write(ref)
    file_ref.close()

    # Save difference file
    file_diff = cdms2.open(parameter.case_id + '/' +
                           parameter.output_file + '_diff.nc', 'w+')
    diff.id = parameter.variables + '(test - reference)'
    file_diff.write(diff)
    file_diff.close()


parser = ACMEParser()
original_parameter = parser.get_parameter(default_vars=False)
parameters = make_parameters(original_parameter)
viewer = OutputViewer(path=parameters[0].case_id, index_name='index name')
add_page_and_top_row(viewer, parameters)

for parameter in parameters:
    viewer.add_group(parameter.case_id)

    reference_data_path = parameter.reference_data_path
    test_data_path = parameter.test_data_path

    var = parameter.variables
    seasons = parameter.season
    ref_name = parameter.ref_name
    test_name = parameter.test_name
    regions = parameter.region

    for season in seasons:
        if hasattr(parameter, 'test_path'):
            filename1 = parameter.test_path
            print filename1
        else:
            test_files = glob.glob(os.path.join(
                test_data_path, '*' + test_name + '*.nc'))
            for filename in fnmatch.filter(test_files, '*' + season + '*'):
                print filename
                filename1 = filename

        if hasattr(parameter, 'reference_path'):
            filename2 = parameter.reference_path
            print filename2
        else:
            ref_files = glob.glob(os.path.join(
                reference_data_path, '*' + ref_name + '*.nc'))
            for filename in fnmatch.filter(ref_files, '*' + ref_name + '_' + season + '*'):
                print filename
                filename2 = filename

        f_mod = cdms2.open(filename1)
        f_obs = cdms2.open(filename2)

        mv1 = acme.process_derived_var(
            var, acme.derived_variables, f_mod, parameter)
        mv2 = acme.process_derived_var(
            var, acme.derived_variables, f_obs, parameter)
        print mv1.units, mv2.units

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

                parameter.output_file = '_'.join(
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
                if hasattr(parameter, 'plot'):
                    parameter.plot(mv2_domain, mv1_domain, diff,
                                   metrics_dict, parameter)
                else:
                    plot(mv2_domain, mv1_domain, diff, metrics_dict, parameter)

                if season is seasons[0]:
                    viewer.add_row('%s %s' % (var, region))
                    viewer.add_col(mv1.long_name)
                viewer.set_row('%s %s' % (var, region))
                files = [parameter.case_id + '/' + parameter.output_file +
                         ext for ext in ['_test.nc', '_ref.nc', '_test.nc']]
                formatted_files = [{'url': f, 'title': f} for f in files]
                viewer.add_col(parameter.case_id + '/' + parameter.output_file +
                               '.png', is_file=True, title=season, other_files=formatted_files)

                save_ncfiles(mv1_domain, mv2_domain, diff, parameter)
                f_mod.close()
                f_obs.close()

        elif mv1.getLevel() and mv2.getLevel():  # for variables with z axis:
            plev = parameter.levels
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

                    parameter.output_file = '_'.join(
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

                    if hasattr(parameter, 'plot'):
                        parameter.plot(mv2_domain, mv1_domain,
                                       diff, metrics_dict, parameter)
                    else:
                        acme_diags.plotting.set5.plot.plot(
                            mv2_domain, mv1_domain, diff, metrics_dict, parameter)

                    if season is seasons[0]:
                        viewer.add_row('%s %s %s' % (
                            var, str(int(plev[ilev])) + 'mb', region))
                        viewer.add_col('Description for %s' % var)
                    viewer.set_row('%s %s %s' %
                                   (var, str(int(plev[ilev])) + 'mb', region))
                    files = [parameter.case_id + '/' + parameter.output_file +
                             ext for ext in ['_test.nc', '_ref.nc', '_test.nc']]
                    formatted_files = [{'url': f, 'title': f} for f in files]
                    viewer.add_col(parameter.case_id + '/' + parameter.output_file +
                                   '.png', is_file=True, title=season, other_files=formatted_files)

                    save_ncfiles(mv1_domain, mv2_domain, diff, parameter)

        else:
            raise RuntimeError(
                "Dimensions of two variables are difference. Abort")

viewer.generate_viewer()
