import os
import copy
import vcs
import cdms2
from genutil import udunits
from acme_diags.driver.utils import get_output_dir


def plot(ref, test, diff, metrics_dict, parameters):
    vcs_canvas = vcs.init()
    ref_test_template = vcs.createtemplate('ref_test_template')
    ref_test_template.blank(["mean", "max", "min", "zvalue", "dataname", "crtime", "ytic2", "xtic2", "xname", "yname", "legend"])

    print test.getAxisList()
    ax = test.getAxis(2)  # Order is T, Y, X
    print dir(ax)
    print ax.getData()
    '''
    ref_test_template.ylabel1.priority = 1
    ref_test_template.ylabel2.priority = 1
    ref_test_template.ymintic1.priority = 1
    ref_test_template.ymintic2.priority = 1
    ref_test_template.yname.priority = 1
    ref_test_template.ytic1.priority = 1
    ref_test_template.ytic2.priority = 1
    ref_test_template.yunits.priority = 1
    ref_test_template.yvalue.priority = 1
    '''

    graph1 = vcs_canvas.createxvsy('ref_test_plot')

    #graph1.datawc_y1 = min(ax.getData())
    #graph1.datawc_y2 = max(ax.getData())
    graph1.datawc_y1 = test.min()
    graph1.datawc_y2 = test.max()
    #graph1.datawc_x1 = 1

    blank_template = vcs.createtemplate()
    blank_template.blank()

    vcs_canvas.plot(test, graph1, ref_test_template)
    vcs_canvas.png('set3vcs.png')


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

def convert_units(var, target_units):
    ''' Converts units of var to target_units.
    var is a cdms.TransientVariable. '''

    if not hasattr(var, 'units') and var.id == 'SST':
        var.units = target_units
    elif not hasattr(var, 'units') and var.id =='ICEFRAC':
        var.units = target_units
        var = 100.0 * var
    elif var.units == 'fraction':
        var = 100.0 * var
        var.units = target_units
    elif var.units == 'mb':
        var.units = target_units
    elif var.units == 'gpm':  # geopotential meter
        var = var / 9.8 / 100  # convert to hecto meter
        var.units = target_units
    else:
        temp = udunits(1.0, var.units)
        coeff, offset = temp.how(target_units)
        var = coeff * var + offset
        var.units = target_units

    return var

if __name__ == '__main__':
    test_pth = '/Users/shaheen2/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'
    ref_pth = '/Users/shaheen2/obs_for_diagnostics/CRU_ANN_climo.nc'

    f = cdms2.open(ref_pth)
    ref = f('TREFHT_LAND')
    f.close()

    f = cdms2.open(test_pth)
    # TREFHT', 'LANDFRAC
    trefht = f('TREFHT')
    landfrac = f('LANDFRAC')
    test = mask_by(convert_units(trefht, target_units="K"), landfrac, low_limit=0.65)
    #test = convert_units(trefht, target_units="K")
    f.close()

    #diff = ref - test
    plot(ref, test, ref, None, None) 
