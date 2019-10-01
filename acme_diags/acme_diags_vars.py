#!/usr/bin/env python
"""
What variables in the passed in file can have E3SM Diagnostics ran on them?
Pass in an E3SM model output file.
It's assumed that this file will have all of the E3SM variables in it.
This is used to get the correct variable names from the derived variables dictionary.
"""
import os
import argparse
import glob
import traceback
import cdms2
import acme_diags
import cdms2.tvariable
from acme_diags.acme_diags_driver import get_parameters
from acme_diags.parser.core_parser import CoreParser
from acme_diags.derivations.acme import derived_variables


# turn off MPI in cdms2 -- not currently supported by e3sm_diags
cdms2.tvariable.HAVE_MPI = False

def main():
    vars_in_e3sm_diags = list_of_vars_in_e3sm_diags()
    vars_with_derived_vars = sorted(check_for_derived_vars(vars_in_e3sm_diags))
    print('Below are the variables needed to run all of the diagnostics in e3sm_diags.')
    print('NOTE: This list doesn\'t include auxiliary variables such as hyam, hybm, PS, etc.')
    print(vars_with_derived_vars)

def list_of_vars_in_user_file():
    """
    Given a path to an nc file, return all of the variables in it.
    """
    #parser = argparse.ArgumentParser()
    #parser.add_argument("path")
    #path = parser.parse_args().path
    # path = DUMMY_FILE_PATH
    path = parser.parse_args().path
    print('Using the file: {}'.format(path))

    if not os.path.exists(path):
        msg = 'The file ({}) does not exist.'.format(path)
        raise RuntimeError(msg)
    with cdms2.open(path) as f:
        return f.variables.keys()

parser = CoreParser()
parser.add_argument('path', default=DUMMY_FILE_PATH, nargs='?')

def list_of_vars_in_e3sm_diags():
    """
    Get a list of all of the variables used in e3sm_diags.
    Open all of the *.cfg files located in acme_diags/acme_diags/driver/default_diags/
    and get all of the 'variables' parameters.
    """

    # Get all of the 'variables' parameter from each file.
    vars_used = []
    #print('hi')
    #args = parser.parse_args()
    #print(dir(args))
    try:
        print('something')
        parameters = get_parameters(parser)
        print('USING USER ARGUMENTS')
    except Exception as e:
        print(e)
        # Looks for these files in their installed location.
        pth = os.path.join(acme_diags.INSTALL_PATH)
        # The first '*' is the folder of the set, the second is the actual file.
        # Ex: {acme_diags.INSTALL_PATH}/lat_lon/lat_lon_model_vs_obs.cfg
        file_paths = [p for p in glob.glob(pth + '*/*.cfg')]
        ## NOT NEEDED:
        ## parser.add_argument('path')  # Needed so the filename can be passed in.
        #parser.add_args_and_values([DUMMY_FILE_PATH])
        parameters = parser.get_other_parameters(files_to_open=file_paths, check_values=False)

    for p in parameters:
        print('p.variables', p.variables)
        vars_used.extend(p.variables)

    # print('vars_used', sorted(list(set(vars_used))))
    # We convert to a set because we only want one of each variable.
    return set(vars_used)

def check_for_derived_vars(e3sm_vars):
    """
    For any of the e3sm_vars which are derived variables, we need
    to check whether any of the original variables are actually in the user's file.

    Ex:
    'PRECT' is a variable in e3sm_vars.
    But it maps to both ('pr',) and ('PRECC', 'PRECL').
    Which one do we use?

    Given a path to a file, we get the vars in that file and
    decided whether to use ('pr',) or ('PRECC', 'PRECL').
    """
    vars_used = []
    vars_in_user_file = set(list_of_vars_in_user_file())
    for var in e3sm_vars:
        if var in derived_variables:
            # Ex: {('PRECC', 'PRECL'): func, ('pr',): func1, ...}.
            vars_to_func_dict = derived_variables[var]
            # Ex: [('pr',), ('PRECC', 'PRECL')].
            possible_vars = vars_to_func_dict.keys()

            var_added = False
            for list_of_vars in possible_vars:
                if not var_added and vars_in_user_file.issuperset(list_of_vars):
                    # All of the variables (list_of_vars) are in the input file.
                    # These are needed.
                    vars_used.extend(list_of_vars)
                    var_added = True
            # If none of the original vars are in the file, just keep this var.
            # This means that it isn't a derived variable in E3SM.
            if not var_added:
                vars_used.append(var)

        else:
            # This var is not a derived variable, it's okay.
            vars_used.append(var)

    return list(set(vars_used))


if __name__ == '__main__':
    main()
