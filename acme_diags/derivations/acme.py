import copy
from numbers import Number
import cdms2
from genutil import udunits


def process_derived_var(var_key, derived_vars_dict, nc_file, parameter):
    ''' Given a key (var_key) to the derived_vars_dict dict, compute and return
     whatever is described in derived_vars_dict[var_key] for the nc_file'''
    if hasattr(parameter, 'derived_variables'):
        _add_user_derived_vars(derived_vars_dict, parameter)
    if var_key in derived_vars_dict:
        return _compute_derived_var(var_key, derived_vars_dict, nc_file)
    else:
        raise RuntimeError(
            'The variable %s was not in the derived variables dictionary' % var_key)


def _add_user_derived_vars(derived_vars_dict, parameter):
    ''' Append parameter.derived_variables to the correct part of derived_vars_dict '''
    for key, user_derived_vars in parameter.derived_variables.iteritems():
        # append the user-defined vars to the already defined ones
        # add to an existing entry, otherwise create a new one
        if key in derived_vars_dict:
            for v in user_derived_vars:
                derived_vars_dict[key][v] = user_derived_vars[v]
        else:
            derived_vars_dict[key] = user_derived_vars


def _compute_derived_var(var_key, derived_vars_dict, nc_file):
    ''' Call the first valid derivation from the derived_vars_dict dict. '''
    derived_vars = derived_vars_dict[var_key]
    # store a list of all inputs visited, so if Exception, we get a good msg.
    derived_var_inputs = []

    # get the first function and inputs from the derived_vars_dict dict
    for inputs, func in derived_vars.iteritems():
        derived_var_inputs.append(inputs)
        # tuples with a single string [ex: ('pr')] become just a string ['pr']
        # are all of the variables (inputs) in the nc_file?
        if isinstance(inputs, str) and inputs in nc_file.variables.keys():
            args = [nc_file(inputs)(squeeze=1)]
            return func(*args)

        elif isinstance(inputs, tuple) and set(inputs).issubset(nc_file.variables.keys()):
            args = [nc_file(var)(squeeze=1) for var in inputs]
            return func(*args)

    # When nc_file is obs, there is var_key in nc_file, i.e. nc_file(var_key)
    if var_key in nc_file.variables.keys():
        return nc_file(var_key)(squeeze=1)

    raise RuntimeError('None of the variables (%s) are in the file: %s' % (
        derived_var_inputs, nc_file.id))


def rename(new_name):
    ''' Given the new name, just return it. '''
    return new_name


def aplusb(var1, var2, target_units=None):
    ''' Returns var1 + var2. If both of their units are not the same,
    it tries to convert both of their units to target_units '''

    if target_units is not None:
        var1 = convert_units(var1, target_units)
        var2 = convert_units(var2, target_units)

    return var1 + var2


def convert_units(var, target_units):
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


def qflxconvert_units(var):
    print var.units
    if var.units == 'kg/m2/s':
        # need to find a solution for units not included in udunits
        # var = convert_units( var, 'kg/m2/s' )
        var = var * 3600.0 * 24  # convert to mm/day
        var.units = 'mm/day'
    return var

def prect(precc, precl):
    """Total precipitation flux = convective + large-scal, """
    var = precc + precl
    var = convert_units(var, "mm/day")
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

def swcfsrf(fsns, fsnsc):
    """Surface shortwave cloud forcing """
    var = fsns - fsnsc
    var.long_name = "Surface shortwave cloud forcing"
    return var

def lwcfsrf(flnsc, flns):
    """Surface longwave cloud forcing """
    var = flnsc - flns
    var.long_name = "Surface longwave cloud forcing"
    return var

def swcf(fsntoa, fsntoac):
    """TOA shortwave cloud forcing """
    var = fsntoa - fsntoac
    var.long_name = "TOA shortwave cloud forcing"
    return var

def lwcf(flntoa, flntoac):
    """TOA longwave cloud forcing """
    var = flntoa - flntoac
    var.long_name = "TOA longwave cloud forcing"
    return var

def fldsc(ts, flnsc):
    """Clearsky Surf LW downwelling flux"""
    var = 5.67e-8 * ts**4 - flnsc
    var.units = "W/m2"
    var.long_name = "Clearsky Surf LW downwelling flux"
    return var

def restom(fsnt, flnt):
    """TOM(top of model) Radiative flux"""
    var = fsnt - flnt
    var.long_name = "TOM(top of model) Radiative flux"
    return var

def restoa(fsnt, flnt):
    """TOA(top of atmosphere) Radiative flux"""
    var = fsnt - flnt
    var.long_name = "TOA(top of atmosphere) Radiative flux"
    return var

# derived_variables is a dictionary to accomodate user specified derived variables. 
# The driver search for available variable keys and functions to calculate derived variable.
# For example 
# In derived_variable there is an entry for 'PRECT':
# If 'PRECT' is not available, but both 'PRECC' and 'PRECT' are available in the netcdf variable keys,
# PRECT is calculated using fuction prect() with precc and precl as inputs.

derived_variables = {
    'PRECT': {
        ('pr'): rename,
        ('PRECC', 'PRECL'): lambda precc, precl: prect(precc, precl)
    },
    'SST': {
        ('SST'): lambda sst: convert_units(sst, target_units="degC"),
        ('TS', 'OCNFRAC'): lambda ts, ocnfrac: mask_by(
            convert_units(ts, target_units="degC"), ocnfrac, low_limit=0.9)
    },
    'PREH2O': {
        ('TMQ'): rename
    },
    'ALBEDO': {
        ('ALBEDO'): rename,
        ('SOLIN', 'FSNTOA'): lambda solin, fsntoa: albedo(solin, fsntoa)
    },
    'ALBEDOC': {
        ('ALBEDOC'): rename,
        ('SOLIN', 'FSNTOAC'): lambda solin, fsntoac: albedoc(solin, fsntoac)
    },
    'SWCF': {
        ('SWCF'): rename,
        ('FSNTOA', 'FSNTOAC'): lambda fsntoa, fsntoac: swcf(fsntoa, fsntoac)
    },
    'SWCFSRF': {
        ('SWCFSRF'): rename,
        ('FSNS', 'FSNSC'): lambda fsns, fsnsc: swcfsrf(fsns, fsnsc)
    },
    'LWCF': {
        ('LWCF'): rename,
        ('FLNTOA', 'FLNTOAC'): lambda flntoa, flntoac: lwcf(flntoa, flntoac)
    },
    'LWCFSRF': {
        ('LWCFSRF'): rename,
        ('FLNSC', 'FLNS'): lambda flns, flnsc: lwcfsrf(flnsc, flns)
    },
    'FLNS': {
        ('FLNS'): rename
    },
    'FLNSC': {
        ('FLNSC'): rename
    },
    'FLDS': {
        ('FLDS'): rename
    },
    'FLDSC': {
        ('FLDSC'): rename,
        ('TS', 'FLNSC'): lambda ts, flnsc: fldsc(ts, flnsc)
    },
    'FSNS': {
        ('FSNS'): rename
    },
    'FSNSC': {
        ('FSNSC'): rename
    },
    'FSDS': {
        ('FSDS'): rename
    },
    'FSDSC': {
        ('FSDSC'): rename
    },
    'FLUT': {
        ('FLUT'): rename
    },
    'FLUTC': {
        ('FLUTC'): rename
    },
    'FSNTOA': {
        ('FSNTOA'): rename
    },
    'FSNTOAC': {
        # Note: CERES_EBAF data in amwg obs sets misspells "units" as "lunits"
        ('FSNTOAC'): rename
    },
    'RESTOM': {
        ('RESTOA'): rename,
        ('FSNT', 'FLNT'): lambda fsnt, flnt: restom(fsnt,flnt)
    },
    'RESTOA': {
        ('RESTOA'): rename,
        ('FSNT', 'FLNT'): lambda fsnt, flnt: restoa(fsnt,flnt)
    },
    'TREFHT_LAND': {
        ('TREFHT_LAND'): rename,
        ('TREFHT', 'LANDFRAC'): lambda trefht, landfrac: mask_by(
            convert_units(trefht, target_units="K"), landfrac, low_limit=0.65)
    },
    'PRECT_LAND': {
        ('PRECIP_LAND'): rename,
        # 0.5 just to match amwg
        ('PRECC', 'PRECL', 'LANDFRAC'): lambda precc, precl, landfrac: mask_by(
           prect(precc , precl), landfrac, low_limit=0.5)
    },
    'Z3': {
        ('Z3'): lambda z3: convert_units(z3, target_units="hectometer")
    },
    'PSL': {
        ('PSL'): lambda psl: convert_units(psl, target_units="mbar")
    },
    'T': {
        ('T'): lambda t: convert_units(t, target_units="K")
    },
    'U': {
        ('U'): lambda u: convert_units(u, target_units="m/s")
    },
    'TREFHT': {
        ('TREFHT'): lambda t: convert_units(t, target_units="K")
    },
    'TREFHT': {
        ('TREFHT'): lambda t: convert_units(t, target_units="K")
    },
    'QFLX': {
        ('QFLX'): lambda qflx: qflxconvert_units(qflx)
    },
    'LHFLX': {
        ('LHFLX'): rename
    },
    'SHFLX': {
        ('SHFLX'): rename
    },
    'TGCLDLWP_OCN': {
        ('TGCLDLWP_OCEAN'): lambda x: convert_units(x, target_units='g/m^2'),
        ('TGCLDLWP', 'OCNFRAC'): lambda tgcldlwp, ocnfrac: mask_by(convert_units(tgcldlwp, target_units="g/m^2"), ocnfrac, low_limit=0.65)
    },
    'PRECT_OCN': {
        ('PRECT_OCEAN'): lambda x: convert_units(x, target_units='mm/day'),
        ('PRECC', 'PRECL', 'OCNFRAC'): lambda a, b, ocnfrac: mask_by(aplusb(a, b, target_units="mm/day"), ocnfrac, low_limit=0.65)
    },
    'PREH2O_OCN': {
        ('PREH2O_OCEAN'): lambda x: convert_units(x, target_units='mm'),
        ('TMQ', 'OCNFRAC'): lambda preh2o, ocnfrac: mask_by(preh2o, ocnfrac, low_limit=0.65)
    },
    'CLDHGH': {
        ('CLDHGH'): lambda cldhgh: convert_units(cldhgh, target_units="%")
    },
    'CLDLOW': {
        ('CLDLOW'): lambda cldlow: convert_units(cldlow, target_units="%")
    },
    'CLDMED': {
        ('CLDMED'): lambda cldmed: convert_units(cldmed, target_units="%")
    },
    'CLDTOT': {
        ('CLDTOT'): lambda cldtot: convert_units(cldtot, target_units="%")
    }
}
