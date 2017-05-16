import os
import copy
import cdutil
import MV2
import genutil
import cdms2
from acme_diags.derivations.default_regions import regions_specs

def findfile(path_name, data_name, season):
    """Locate file name based on data_name and season."""
    dir_files = os.listdir(path_name)
    for filename in dir_files:
        if filename.startswith(data_name) and season in filename:
            return path_name+filename
    raise IOError(
        "No file found based on given path_name and data_name")
       
def hybrid_to_plevs(var, hyam, hybm, ps, plev):
    """Convert from hybrid pressure coordinate to desired pressure level(s)."""
    p0 = 1000. # mb
    ps = ps /100. # convert unit from 'Pa' to mb
    levels_orig = cdutil.vertical.reconstructPressureFromHybrid(
        ps, hyam, hybm, p0)
    levels_orig.units = 'mb'
    # Make sure z is positive down
    if var.getLevel()[0] > var.getLevel()[-1]:
        var = var(lev=slice(-1,None,-1)) 
        levels_orig = levels_orig(lev=slice(-1,None,-1)) 
    var_p = cdutil.vertical.logLinearInterpolation(
        var(squeeze=1), levels_orig(squeeze=1), plev)

    return var_p

def pressure_to_plevs(var, plev):
    """Convert from pressure coordinate to desired pressure level(s)."""
    # Construct pressure level for interpolation
    var_plv = var.getLevel()
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
        var = var(lev=slice(-1,None,-1)) 
        levels_orig = levels_orig(lev=slice(-1,None,-1)) 
    var_p = cdutil.vertical.logLinearInterpolation(
        var(squeeze=1), levels_orig(squeeze=1), plev)
    return var_p


def select_region(region, var1, var2, land_frac,ocean_frac,parameter):
    """Select desired regions from transient variables."""
    domain = None
    # if region != 'global':
    if region.find('land') != -1 or region.find('ocean') != -1:
        if region.find('land') != -1:
            land_ocean_frac = land_frac
        elif region.find('ocean') != -1:
            land_ocean_frac = ocean_frac
        region_value = regions_specs[region]['value']
        print 'region_value', region_value
    
        var1_domain = mask_by(
            var1, land_ocean_frac, low_limit=region_value)
        var2_domain = var2.regrid(
            var1.getGrid(), parameter.regrid_tool, parameter.regrid_method)
        var2_domain = mask_by(
            var2_domain, land_ocean_frac, low_limit=region_value)
    else:
        var1_domain = var1
        var2_domain = var2
    
    try:
        # if region.find('global') == -1:
        domain = regions_specs[region]['domain']
        print domain
    except:
        print ("no domain selector")
    var1_domain = var1_domain(domain)
    var2_domain = var2_domain(domain)
    var1_domain.units = var1.units
    var2_domain.units = var1.units
    
    return var1_domain, var2_domain


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
    else:
        mv_grid = mv2.getGrid()
        mv2_reg = mv2
        mv1_reg = mv1.regrid(mv_grid, regridTool=regrid_tool,
                             regridMethod=regrid_method)
    return mv1_reg, mv2_reg

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


def save_ncfiles(set_num, test, ref, diff, parameter):
    """Saves the test, reference, and difference nc files."""
    # Save files being plotted
    # Set cdms preferences - no compression, no shuffling, no complaining
    cdms2.setNetcdfDeflateFlag(1)
    # 1-9, min to max - Comes at heavy IO (read/write time cost)
    cdms2.setNetcdfDeflateLevelFlag(0)
    cdms2.setNetcdfShuffleFlag(0)
    cdms2.setCompressionWarnings(0)  # Turn off warning messages
    # Save test file
    pth = get_output_dir(set_num, parameter)
    file_test = cdms2.open(pth + '/' + parameter.output_file + '_test.nc', 'w+')
    test.id = parameter.var_id
    file_test.write(test)
    file_test.close()

    # Save reference file
    file_ref = cdms2.open(pth + '/' + parameter.output_file + '_ref.nc', 'w+')
    ref.id = parameter.var_id
    file_ref.write(ref)
    file_ref.close()

    # Save difference file
    file_diff = cdms2.open(pth + '/' + parameter.output_file + '_diff.nc', 'w+')
    diff.id = parameter.var_id + '(test - reference)'
    file_diff.write(diff)
    file_diff.close()

def get_output_dir(set_num, parameter):
    """Get the directory of where to save the outputs for a run."""
    pth = os.path.join(parameter.results_dir, 'set{}'.format(set_num), parameter.case_id)
    if not os.path.exists(pth):
        os.makedirs(pth)
    return pth
    
