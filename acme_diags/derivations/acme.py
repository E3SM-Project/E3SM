from numbers import Number
from unidata import udunits

def process_derived_var(var_key, derived_vars_dict, nc_file, region):
    ''' Given a key (var_key) to the derived_vars_dict dict, compute and return
     whatever is described in derived_vars_dict[var_key] for the nc_file'''
    if var_key in derived_vars_dict.keys():
        inputs, func = _get_correct_derivation(var_key, derived_vars_dict, nc_file)
        # get all of the variables from nc_file
        args = [nc_file(var, region)(squeeze=1) for var in inputs]
        #args = [nc_file(var)(squeeze=1) for var in inputs]
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

derived_variables = {
    'PRECT': [
        (['pr'],  rename),
        (['PRECC','PRECL'], lambda a, b: aplusb(a, b, target_units="mm/day"))
    ],
}
