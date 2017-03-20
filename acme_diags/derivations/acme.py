from numbers import Number
from unidata import udunits
import cdms2
import copy


def process_derived_var(var_key, derived_vars_dict, nc_file):
    ''' Given a key (var_key) to the derived_vars_dict dict, compute and return
     whatever is described in derived_vars_dict[var_key] for the nc_file'''
    if var_key in derived_vars_dict.keys():
    #if var_key:
        inputs, func = _get_correct_derivation(var_key, derived_vars_dict, nc_file)
        # get all of the variables from nc_file
        args = [nc_file(var)(squeeze=1) for var in inputs]
        return func(*args)
    else:
        raise RuntimeError('The variable %s was not in the derived variables dictionary' % var_key)

def _get_correct_derivation(var_key, derived_vars_dict, nc_file):
    ''' Get the first valid derivation from the derived_vars_dict dict. '''
    derived_var_list = derived_vars_dict[var_key]
    derived_var_inputs = []  # store a list of all inputs visited, so if Exception, we get a good msg.

    # When nc_file is obs, there is var_key in nc_file, i.e. nc_file(var_key)
    if var_key in nc_file.variables.keys():
        return [var_key], lambda x: x  # var_key needs to be in a list

    # get the first function and inputs from the derived_vars_dict dict
    for inputs, func in derived_var_list:
        derived_var_inputs.append(inputs)
        # are all of the variables (inputs) in the nc_file?
        if set(inputs).issubset(nc_file.variables.keys()):
            return inputs, func

    raise RuntimeError('None of the variables (%s) are in the file: %s' % (derived_var_inputs, nc_file.id))

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

    temp = udunits(1.0, var.units)
    coeff, offset= temp.how(target_units)

    var = coeff*var + offset
    var.units = target_units

    return var

def mask_by( input_var, maskvar, low_limit=None, high_limit=None ):
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

#def qflx_to_lhflx( var ):
#    """computes latent heat flux from q flux ."""
#
#    qflx = _convert_units( var, 'kg/(m^2 s)' )
#    lv = 2.501e6              # latent heat of evaporation units J/kg 
#    lhflx = qflx *lv
#    lhflx.units = "W/m^2" 
#    lhflx.long_name = "Surf latent heat flux"
#    return lhflx

def qflx_convert_units (var):
    print "testtest"
    var = _convert_units( var, 'kg/(m^2 s)' )
    var = var *3600.0*24  #convert to mm/day
    print '%%%%%%%%%%%%%%%%%%%%%%%%%'
    var.units = 'mm/day'
    return var


derived_variables = {
    'PRECT': [
        (['pr'],  rename),
        (['PRECC','PRECL'], lambda a, b: aplusb(a, b, target_units="mm/day"))
    ],
    'SST': [
        (['TS', 'OCNFRAC'], lambda ts, ocnfrac: mask_by(_convert_units(ts, target_units="degC"), ocnfrac, low_limit = 0.9))
    ],
    'PREH2O': [
        (['TMQ'], rename)
    ],
    'TREFHT_LAND':[ 
        (['TREFHT', 'LANDFRAC'], lambda trefht, landfrac: mask_by(_convert_units(trefht, target_units="K"), landfrac, low_limit = 0.65)),
        (['TREFHT_LAND'],rename)
    ],
    'PRECT_LAND':[ 
        (['PRECC','PRECL', 'LANDFRAC'], lambda a, b , landfrac: mask_by(aplusb(a, b, target_units="mm/day"), landfrac, low_limit = 0.65) ),
        (['PRECIP_LAND'],rename)
    ],
    'Z3': [
        (['Z3'], lambda z3: _convert_units(z3, target_units="hectometer"))
    ],
    'PSL': [
        (['PSL'], lambda psl: _convert_units(psl, target_units="hectopascal"))
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
    'QFLX': [
        (['QFLX'], lambda qflx: qflx_convert_units(qflx))
    ],
    'LHFLX': [
        (['LHFLX'], lambda lhflx: _convert_units(lhflx, target_units="W/m^2")),
    ]
}

