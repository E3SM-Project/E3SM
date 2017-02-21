from numbers import Number
from unidata import udunits

def process_derived_var(var_key, nc_file):
    ''' Given a key to the derived_variables dict, compute and return
     whatever is described in derived_variables[var_key] '''
    # This function will interface with derived_variables and any other dict
    # in the derivations folder.
    if var_key in derived_variables.keys():
        inputs, func = _get_correct_derivation(var_key, nc_file)
        # get all of the variables from nc_file
        args = [nc_file(var)(squeeze=1) for var in inputs]
        return func(*args)
    else:
        raise RuntimeError('The variable %s was not in the dictionary' % var)

def _get_correct_derivation(var_key, nc_file):
    ''' Get the first valid derivation from the derived_variables dict. '''
    derived_var_list = derived_variables[var_key]
    # get the first function and args from
    for inputs, func in derived_var_list:
        # are all of the variables (inputs) in the nc_file?
        if set(inputs).issubset(nc_file.variables.keys()):
            return inputs, func

def rename(new_name):
    ''' Given the new name, just return it. '''
    return new_name

def aplusb(var1, var2, target_units=None):
    ''' Returns var1 + var2. If both of their units are not the same,
    it tries to convert both of their units to target_units '''

    if target_units is not None:
        var1, var2 = _convert_units(var1, var2, target_units)

    return var1 + var2

def _convert_units(var1, var2, target_units):
    ''' Converts units of var1 and var2 to target_units.
    var1 and var2 are cdms.TransientVariable. '''
    #return udunits(target_units, var1.units)
    # temp is of type udunits
    temp = udunits(1.0, var1.units)
    coeff, offset= temp.how(target_units)

    var1 = coeff*var1 + offset
    var2 = coeff*var2 + offset

    var1.units = target_units
    var2.units = target_units

    return var1, var2

derived_variables = {
    'PRECT': [
        (['pr'],  rename),
        (['PRECC','PRECL'], lambda a, b: aplusb(a, b, target_units="mm/day"))
    ],
}
