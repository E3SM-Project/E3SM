def get_var_from_file_by_derived_var(nc_file, var_key):
    ''' Given a key to the derived_variables dict, compute and return
     whatever is described in derived_variables[var_key] '''
    # This function will interface with derived_variables and any other dict
    # in the derivations folder.
    if var_key in derived_variables.keys():
        derived_var_list = derived_variables[var_key]
        # What to do with multiple inputs???
        for inputs, func in derived_var_list:
            args = _get_cdms_var_from_nc_file(inputs, nc_file)
            if args != []:
                return func(*args)
    else:
        raise RuntimeError('The variable %s was not in the dictionary' % var)

def _get_cdms_var_from_nc_file(inputs, nc_file):
    ''' Given a list of string variable inputs, get
    the actual cdms variables from the nc_file'''
    args = []
    for var in inputs:
        if var in nc_file.variables.keys():
            arg += [nc_file(var)(squeeze=1)]
    return args

def rename(new_name):
    ''' Given the new name, return it. '''
    return new_name


derived_variables = {
    'PRECT': [
        (['pr'],  rename),
        (['PRECC','PRECL'], aplusb(a,b,units="mm/day"))
    ],
}
