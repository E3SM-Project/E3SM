from __future__ import print_function

import os
import copy
import pwd
import grp
import cdutil
import MV2
import genutil
import cdms2
from acme_diags import container
from acme_diags.derivations.default_regions import regions_specs

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))

def monotonically_decreasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))

def monotonically_increasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))

def monotonic(L):
        return monotonically_increasing(L) or monotonically_decreasing(L)


def adjust_time_from_time_bounds(var):
    """
    Redefine time to be in the middle of the time interval, and rewrite 
    the time axis. This is important for data where the absolute time doesn't fall in the middle of the time interval, such as E3SM, the time was recorded at the end of each time Bounds.
    """
    var_time = var.getTime()
    tbounds = var_time.getBounds()
    var_time[:] = 0.5*(tbounds[:,0]+tbounds[:,1])
    var_time_absolute = var_time.asComponentTime()
    time2 = cdms2.createAxis(var_time)
    time2.designateTime()
    #.designateTime() needs to be set before attributes changes.
    time2.units = var_time.units
    time2.calendar = var_time.calendar
    time2.setBounds(tbounds)
    #time2.calendar = cdtime.NoLeapCalendar
    time2.id = 'time'
    var.setAxis(0,time2)
#    cdutil.setTimeBoundsMonthly(var)
 
    return var


def get_name_and_yrs(parameters, dataset, season=''):
    """
    Given either test or ref data, get the name of the data
    (test_name or reference_name), along with the years averaged.
    """
    if dataset.test:
        if parameters.short_test_name:
            name_yrs = parameters.short_test_name
        else:
            name_yrs = parameters.test_name
    else:
        if parameters.short_ref_name:
            name_yrs = parameters.short_ref_name
        else:
            # TODO: Or is this ref_name?
            name_yrs = parameters.reference_name

    if dataset.is_climo():
        try:
            yrs_averaged = dataset.get_attr_from_climo('yrs_averaged', season)
            name_yrs = '{} ({})'.format(name_yrs, yrs_averaged)
        except:
            print("No 'yrs_averaged' exists in the global attributes.")
    else:
        start_yr, end_yr = dataset.get_start_and_end_years()
        yrs_averaged = '{}-{}'.format(start_yr, end_yr)
        name_yrs = '{} ({})'.format(name_yrs, yrs_averaged)
    
    return name_yrs


def convert_to_pressure_levels(mv, plevs, dataset, var, season):
    """
    Given either test or reference data with a z-axis,
    convert to the desired pressure levels.
    """
    mv_plv = mv.getLevel()
    # var(time,lev,lon,lat) convert from hybrid level to pressure
    if mv_plv.long_name.lower().find('hybrid') != -1:
        extra_vars = ['hyam', 'hybm', 'PS']
        hyam, hybm, ps = dataset.get_extra_variables_only(var, season, extra_vars=extra_vars)
        mv_p = hybrid_to_plevs(mv, hyam, hybm, ps, plevs)

    # levels are pressure levels
    elif mv_plv.long_name.lower().find('pressure') != -1 or mv_plv.long_name.lower().find('isobaric') != -1:
        mv_p = pressure_to_plevs(mv, plevs)

    else:
        raise RuntimeError(
            "Vertical level is neither hybrid nor pressure. Aborting.")
            
    return mv_p


def hybrid_to_plevs(var, hyam, hybm, ps, plev):
    """Convert from hybrid pressure coordinate to desired pressure level(s)."""
    p0 = 1000.  # mb
    ps = ps / 100.  # convert unit from 'Pa' to mb
    levels_orig = cdutil.vertical.reconstructPressureFromHybrid(
        ps, hyam, hybm, p0)
    levels_orig.units = 'mb'
    # Make sure z is positive down
    if var.getLevel()[0] > var.getLevel()[-1]:
        var = var(lev=slice(-1, None, -1))
        levels_orig = levels_orig(lev=slice(-1, None, -1))
    var_p = cdutil.vertical.logLinearInterpolation(
        var(squeeze=1), levels_orig(squeeze=1), plev)

    return var_p


def pressure_to_plevs(var, plev):
    """Convert from pressure coordinate to desired pressure level(s)."""
    # Construct pressure level for interpolation
    var_plv = var.getLevel()
    if var_plv.units == 'Pa':
        var_plv[:] = var_plv[:]/100.0 #convert Pa to mb
    levels_orig = MV2.array(var_plv[:])
    levels_orig.setAxis(0, var_plv)
    # grow 1d levels_orig to mv dimention
    var, levels_orig = genutil.grower(var, levels_orig)
    # levels_orig.info()
    # logLinearInterpolation only takes positive down plevel:
    # "I :      interpolation field (usually Pressure or depth)
    # from TOP (level 0) to BOTTOM (last level), i.e P value
    # going up with each level"
    if var.getLevel()[0] > var.getLevel()[-1]:
        var = var(lev=slice(-1, None, -1))
        levels_orig = levels_orig(lev=slice(-1, None, -1))
    var_p = cdutil.vertical.logLinearInterpolation(
        var(squeeze=1), levels_orig(squeeze=1), plev)

    return var_p


def select_region(region, var, land_frac, ocean_frac, parameter):
    """Select desired regions from transient variables."""
    domain = None
    # if region != 'global':
    if region.find('land') != -1 or region.find('ocean') != -1:
        if region.find('land') != -1:
            land_ocean_frac = land_frac
        elif region.find('ocean') != -1:
            land_ocean_frac = ocean_frac
        region_value = regions_specs[region]['value']

        land_ocean_frac = land_ocean_frac.regrid(
            var.getGrid(), regridTool=parameter.regrid_tool, regridMethod=parameter.regrid_method)

        var_domain = mask_by(
            var, land_ocean_frac, low_limit=region_value)
    else:
        var_domain = var

    try:
        # if region.find('global') == -1:
        domain = regions_specs[region]['domain']
        #print('Domain: ', domain)
    except:
        pass
        #print("No domain selector.")

    var_domain_selected = var_domain(domain)
    var_domain_selected.units = var.units

    return var_domain_selected


def regrid_to_lower_res(mv1, mv2, regrid_tool, regrid_method):
    """Regrid transient variable toward lower resolution of two variables."""

    axes1 = mv1.getAxisList()
    axes2 = mv2.getAxisList()

    # use nlat to decide data resolution, higher number means higher data
    # resolution. For the difference plot, regrid toward lower resolution
    if len(axes1[1]) <= len(axes2[1]):
        mv_grid = mv1.getGrid()
        mv1_reg = mv1
        mv2_reg = mv2.regrid(mv_grid, regridTool=regrid_tool,
                             regridMethod=regrid_method)
        mv2_reg.units = mv2.units

    else:
        mv_grid = mv2.getGrid()
        mv2_reg = mv2
        mv1_reg = mv1.regrid(mv_grid, regridTool=regrid_tool,
                             regridMethod=regrid_method)
        mv1_reg.units = mv1.units

    return mv1_reg, mv2_reg


def mask_by(input_var, maskvar, low_limit=None, high_limit=None):
    """masks a variable var to be missing except where maskvar>=low_limit and maskvar<=high_limit.
    None means to omit the constrint, i.e. low_limit = -infinity or high_limit = infinity.
    var is changed and returned; we don't make a new variable.
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


def save_ncfiles(set_num, test, ref, diff, parameter):
    """
    Saves the test, reference, and difference
    data being plotted as nc files.
    """
    if parameter.save_netcdf:
        # Save files being plotted
        # Set cdms preferences - no compression, no shuffling, no complaining
        cdms2.setNetcdfDeflateFlag(1)
        # 1-9, min to max - Comes at heavy IO (read/write time cost)
        cdms2.setNetcdfDeflateLevelFlag(0)
        cdms2.setNetcdfShuffleFlag(0)
        cdms2.setCompressionWarnings(0)  # Turn off warning messages

        pth = get_output_dir(set_num, parameter)

        # Save test file
        test.id = parameter.var_id
        test_pth = os.path.join(pth, parameter.output_file + '_test.nc')
        with cdms2.open(test_pth, 'w+') as file_test:
            file_test.write(test)

        # Save reference file
        ref.id = parameter.var_id
        ref_pth = os.path.join(pth, parameter.output_file + '_ref.nc')
        with cdms2.open(ref_pth, 'w+') as file_ref:
            file_ref.write(ref)

        # Save difference file
        diff.id = parameter.var_id + '(test - reference)'
        diff_pth = os.path.join(pth, parameter.output_file + '_diff.nc')
        with cdms2.open(diff_pth, 'w+') as file_diff:
            file_diff.write(diff)


def get_output_dir(set_num, parameter, ignore_container=False):
    """
    Get the directory of where to save the outputs for a run.
    If ignore_container is True and the software is being ran in a container,
      get the path that the user passed in.
    """
    results_dir = parameter.results_dir
    if ignore_container and container.is_container():
        results_dir = parameter.orig_results_dir
    
    pth = os.path.join(results_dir,
                       '{}'.format(set_num), parameter.case_id)
    if not os.path.exists(pth):
        # When running diags in parallel, sometimes another process will create the dir.
        try:
            os.makedirs(pth, 0o755)
        except OSError as e:
            if e.errno != os.errno.EEXIST:
                raise

    return pth
