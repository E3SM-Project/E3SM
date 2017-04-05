import copy
from numbers import Number
import cdms2
from genutil import udunits


def process_derived_var(var_key, derived_vars_dict, nc_file, parameter):
    ''' Given a key (var_key) to the derived_vars_dict dict, compute and return
     whatever is described in derived_vars_dict[var_key] for the nc_file'''
    if hasattr(parameter, 'derived_variables'):
        _add_user_derived_vars(derived_vars_dict, parameter)
    if var_key in derived_vars_dict.keys():
        inputs, func = _get_correct_derivation(
            var_key, derived_vars_dict, nc_file)
        # get all of the variables from nc_file
        args = [nc_file(var)(squeeze=1) for var in inputs]
        return func(*args)
    else:
        raise RuntimeError(
            'The variable %s was not in the derived variables dictionary' % var_key)


def _add_user_derived_vars(derived_vars_dict, parameter):
    for k, v in parameter.derived_variables.iteritems():
        derived_vars_dict[k] = v


def _get_correct_derivation(var_key, derived_vars_dict, nc_file):
    ''' Get the first valid derivation from the derived_vars_dict dict. '''
    derived_var_list = derived_vars_dict[var_key]
    # store a list of all inputs visited, so if Exception, we get a good msg.
    derived_var_inputs = []

    # get the first function and inputs from the derived_vars_dict dict
    for inputs, func in derived_var_list:
        derived_var_inputs.append(inputs)
        # are all of the variables (inputs) in the nc_file?
        if set(inputs).issubset(nc_file.variables.keys()):
            return inputs, func

    # When nc_file is obs, there is var_key in nc_file, i.e. nc_file(var_key)
    if var_key in nc_file.variables.keys():
        return [var_key], lambda x: x  # var_key needs to be in a list

    raise RuntimeError('None of the variables (%s) are in the file: %s' % (
        derived_var_inputs, nc_file.id))


def rename(new_name):
    ''' Given the new name, just return it. '''
    return new_name


def aplusb(var1, var2, target_units=None):
    ''' Returns var1 + var2. If both of their units are not the same,
    it tries to convert both of their units to target_units '''

    if target_units is not None:
        var1 = _convert_units(var1, target_units)
        var2 = _convert_units(var2, target_units)

    return var1 + var2


def _convert_units(var, target_units):
    ''' Converts units of var to target_units.
    var is a cdms.TransientVariable. '''

    if not hasattr(var, 'units') and var.id == 'SST':
        var.units = target_units
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


def qflx_convert_units(var):
    print var.units
    if var.units == 'kg/m2/s':
        # need to find a solution for units not included in udunits
        #var = _convert_units( var, 'kg/m2/s' )
        var = var * 3600.0 * 24  # convert to mm/day
        var.units = 'mm/day'
    return var

def prect(precc, precl):
    """Total precipitation flux = convective + large-scal, """
    var = precc + precl
    var = _convert_units(var, "mm/day") 
    var.long_name = "Total precipitation rate (convective + large-scale)"
    return var
    
def albedo(solin, fsntoa):
    """TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension"""
    var = (solin - fsntoa) / solin
    var.units = "dimensionless"
    var.long_name = "TOA albedo"
    return var

def albedoc(solin, fsntoa):
    """TOA (top-of-atmosphere) albedo clear-sky, (solin - fsntoac) / solin, unit is nondimension"""
    var = (solin - fsntoa) / solin
    var.units = "dimensionless"
    var.long_name = "TOA albedo clear-sky"
    return var

# derived_variables is a list
derived_variables = {
    'PRECT': [
        (['pr'], rename),
        (['PRECC', 'PRECL'], lambda precc, precl: prect(precc , precl))
    ],
    'SST': [
        (['SST'], lambda sst: _convert_units(sst, target_units="degC")),
        (['TS', 'OCNFRAC'], lambda ts, ocnfrac: mask_by(
            _convert_units(ts, target_units="degC"), ocnfrac, low_limit=0.9))
    ],
    'PREH2O': [
        (['TMQ'], rename)
    ],
    'ALBEDO': [
        (['ALBEDO'], rename),
        (['SOLIN', 'FSNTOA'], lambda solin, fsntoa: albedo(solin, fsntoa))
    ],
    'ALBEDOC': [
        (['ALBEDOC'], rename),
        (['SOLIN', 'FSNTOAC'], lambda solin, fsntoac: albedoc(solin, fsntoac))
    ],
    'SWCF': [
        (['SWCF'], rename),
        (['FSNTOA', 'FSNTOAC'], lambda fsntoa, fsntoac: fsntoa - fsntoac)
    ],
    'SWCFSRF': [
        (['SWCFSRF'], rename),
        (['FSNS', 'FSNSC'], lambda fsns, fsnsc: fsns - fsnsc)
    ],
    'LWCF': [
        (['LWCF'], rename),
        (['FLNTOA', 'FLNTOAC'], lambda flntoa, flntoac: flntoa - flntoac)
    ],
    'LWCFSRF': [
        (['LWCFSRF'], rename),
        (['FLNSC', 'FLNS'], lambda flns, flnsc: flnsc - flns)
    ],
    'FLNS': [
        (['FLNS'], rename)
    ],
    'FLNSC': [
        (['FLNSC'], rename)
    ],
    'FLDS': [
        (['FLDS'], rename)
    ],
    'FLDSC': [
        (['FLDSC'], rename),
        (['TS', 'FLNSC'], lambda ts, flnsc: 5.67e-8 * ts**4 - flnsc)
    ],
    'FSNS': [
        (['FSNS'], rename)
    ],
    'FSNSC': [
        (['FSNSC'], rename)
    ],
    'FSDS': [
        (['FSDS'], rename)
    ],
    'FSDSC': [
        (['FSDSC'], rename)
    ],
    'FLUT': [
        (['FLUT'], rename)
    ],
    'FLUTC': [
        (['FLUTC'], rename)
    ],
    'FSNTOA': [
        (['FSNTOA'], rename)
    ],
    'FSNTOAC': [
        (['FSNTOAC'], rename)
    ],
    'RESTOM': [
        (['RESTOA'], rename),
        (['FSNT', 'FLNT'], lambda fsnt, flnt: fsnt - flnt)
    ],
    'RESTOA': [
        (['RESTOA'], rename),
        (['FSNTOA', 'FLUT'], lambda fsntoa, flut: fsntoa - flut)
    ],
    'TREFHT_LAND': [
        (['TREFHT', 'LANDFRAC'], lambda trefht, landfrac: mask_by(
            _convert_units(trefht, target_units="K"), landfrac, low_limit=0.65)),
        (['TREFHT_LAND'], rename)
    ],
    'PRECT_LAND': [
        (['PRECC', 'PRECL', 'LANDFRAC'], lambda a, b, landfrac: mask_by(aplusb(
            a, b, target_units="mm/day"), landfrac, low_limit=0.5)),  # 0.5 just to match amwg
        (['PRECIP_LAND'], rename)
    ],
    'Z3': [
        (['Z3'], lambda z3: _convert_units(z3, target_units="hectometer"))
    ],
    'PSL': [
        (['PSL'], lambda psl: _convert_units(psl, target_units="mbar"))
    ],
    'T': [
        (['T'], lambda t: _convert_units(t, target_units="K"))
    ],
    'U': [
        (['U'], lambda u: _convert_units(u, target_units="m/s"))
    ],
    'TREFHT': [
        (['TREFHT'], lambda t: _convert_units(t, target_units="K"))
    ],
    'TREFHT': [
        (['TREFHT'], lambda t: _convert_units(t, target_units="K"))
    ],
    'QFLX': [
        (['QFLX'], lambda qflx: qflx_convert_units(qflx))
    ],

    'LHFLX': [
        (['LHFLX'], rename)
    ],
    'SHFLX': [
        (['SHFLX'], rename)
    ],
    'TGCLDLWP_OCN': [(['TGCLDLWP_OCEAN'], (lambda x: _convert_units(x, target_units='g/m^2'))),
                     (['TGCLDLWP', 'OCNFRAC'], lambda tgcldlwp, ocnfrac: mask_by(_convert_units(tgcldlwp, target_units="g/m^2"), ocnfrac, low_limit=0.65))],
    'PRECT_OCN': [(['PRECT_OCEAN'], (lambda x: _convert_units(x, target_units='mm/day'))),
                  (['PRECC', 'PRECL', 'OCNFRAC'], lambda a, b, ocnfrac: mask_by(aplusb(a, b, target_units="mm/day"), ocnfrac, low_limit=0.65))],
    'PREH2O_OCN': [(['PREH2O_OCEAN'], (lambda x: _convert_units(x, target_units='mm'))),
                   (['TMQ', 'OCNFRAC'], lambda preh2o, ocnfrac: mask_by(preh2o, ocnfrac, low_limit=0.65))],
    'CLDHGH': [
        (['CLDHGH'], lambda cldhgh: _convert_units(cldhgh, target_units="%")),
    ],
    'CLDLOW': [
        (['CLDLOW'], lambda cldlow: _convert_units(cldlow, target_units="%")),
    ],
    'CLDMED': [
        (['CLDMED'], lambda cldmed: _convert_units(cldmed, target_units="%")),
    ],
    'CLDTOT': [
        (['CLDTOT'], lambda cldtot: _convert_units(cldtot, target_units="%")),
    ],
#    'CLDHGH_VISIR': [
#        (['CLDHGH'], rename),
#        (['CLDHGH_VISIR'], rename),
#    ],
#    'CLDLOW_VISIR': [
#        (['CLDLOW'], rename),
#        (['CLDLOW_VISIR'], rename),
#    ],
#    'CLDMED_VISIR': [
#        (['CLDMED'], rename),
#        (['CLDMED_VISIR'], rename),
#    ],
#    'CLDTOT_VISIR': [
#        (['CLDTOT'], rename),
#        (['CLDTOT_VISIR'], rename),
#    ],
}
