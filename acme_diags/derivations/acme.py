def get_something(var_key):
    ''' Given a key to the derived_variables dict, compute and return
     whatever is described in derived_variables[var_key]'''
    # This function will interface with derived_variables and any other dict
    # in the derivations folder.
    if var_key in derived_variables.keys():
        derived_var_list = derived_variables[var_key]
        # What to do with multiple inputs???
        for inputs, func in derived_var_list:
            return func(*inputs)
    else:
        raise RuntimeError('The variable %s was not in the dictionary' % var)

def rename(new_name):
    ''' Given the new name, return it. '''
    return new_name

derived_variables = {
    'PRECT': [
        (['pr'],  lambda x: rename(x)),
        (['pr2'],  lambda x: rename(x)),
        #(['PRECC','PRECL'], (lambda a,b,units="mm/day": aplusb(a,b,units)) ))
    ],
}


print get_something('PRECT')
