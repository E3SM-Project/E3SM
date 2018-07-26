"""
Get a variable from input data (either reference or test data).
This data can either be climatology files or timeseries files.
Derived variables are also supported.
"""
import os
import collections
import traceback
import cdms2
import acme_diags.derivations.acme
from acme_diags.driver import utils


class Dataset():
    def __init__(self, parameters, ref=False, test=False, derived_vars={}):
        self.parameters = parameters
        self.ref = ref
        self.test = test
        self.derived_vars = derived_vars

        if self.ref is False and self.test is False:
            msg = 'Both ref and test cannot be False. One must be True.'
            raise RuntimeError(msg)
        elif self.ref is True and self.test is True:
            msg = 'Both ref and test cannot be True. Only one must be True.'
            raise RuntimeError(msg)
        
        if not self.derived_vars:
            # Use the default derived variables.
            self.derived_vars = acme_diags.derivations.acme.derived_variables

        if hasattr(self.parameters, 'derived_variables'):
            self._add_user_derived_vars()


    def _add_user_derived_vars(self):
        """
        If the user-defined derived variables is in the input parameters, append
        parameters.derived_variables to the correct part of self.derived_vars.
        """
        key_val_pairs = self.parameters.derived_variables.items()
        for derived_var, original_vars in list(key_val_pairs):
            # Append the user-defined vars to the already defined ones.
            if derived_var in self.derived_vars:
                # Put user-defined derived vars first in the OrderedDict.
                # That's why we create a new one.
                new_dict = collections.OrderedDict(original_vars)
                # Add all of the default derived vars to the end of new_dict.
                for k in self.derived_vars[derived_var]:
                    # Don't overwrite the user-defined var with a default derived var.
                    if k in new_dict:
                        continue
                    new_dict[k] = self.derived_vars[derived_var][k]
                self.derived_vars[derived_var] = new_dict
            # Otherwise, this is a new derived var, so add it as a new entry.
            else:
                self.derived_vars[derived_var] = original_vars


    def get_variable(self, var, season, extra_vars=[]):
        """
        For a given season, get the variable and any extra variables.
        These variables can either be from the test data or reference data.
        """
        self.var = var
        self.season = season
        self.extra_vars = extra_vars

        if not self.var:
            raise RuntimeError('Variable is invalid.')
        if not self.season:
            raise RuntimeError('Season is invalid.')

        # We need to make two decisions:
        # 1) Are the files being used reference or test data?
        #    - This is done with self.ref and self.test.
        # 2) Are the files being used climo or timeseries files?
        #   - This is done with the ref_timeseries_input and test_timeseries_input parameters.
        if self.ref and getattr(self.parameters, 'ref_timeseries_input', False) is True:
            # Get the reference variable from timeseries files.
            data_path = self.parameters.ref_data_path
            start_yr = self.parameters.parameters.ref_start_yr
            end_yr = self.parameters.ref_end_yr
            return self._get_timeseries_var(data_path, start_yr, end_yr)
            
        elif self.test and getattr(self.parameters, 'test_timeseries_input', False) is True:
            # Get the test variable from timeseries files.
            data_path = self.parameters.test_data_path
            start_yr = self.parameters.test_start_yr
            end_yr = self.parameters.test_end_yr
            return self._get_timeseries_var(data_path, start_yr, end_yr)

        elif self.ref:
            # Get the reference variable from climo files.
            filename = utils.get_ref_filename(self.parameters, self.season)
            return self._get_climo_var(filename)

        elif self.test:
            # Get the test variable from climo files.
            filename = utils.get_test_filename(self.parameters, self.season)
            return self._get_climo_var(filename)

        else:
            msg = '''
            Error when determining what kind (ref or test) of variable to get and
            where to get it from (climo or timeseries files).
            '''
            raise RuntimeError(msg)
    

    def get_extra_variables_only(self, var, season, extra_vars):
        """
        For a given season, get only the extra variables.
        These can either be from the test data or reference data.
        """
        # Do something like this:
        # stuff = self.get_variable(var, season, extra_vars)
        # return stuff[1:]
        pass


    def _get_climo_var(self, filename):
        """
        For a given season (self.season) and climo
        input data, get the variable (self.var).

        If self.extra_vars is also defined, get them as well.
        """
        vars_to_get = [self.var]
        vars_to_get.extend(self.extra_vars)
        variables = []

        with cdms2.open(filename) as data_file:
            for var in vars_to_get:
                # If it's a derived var, get that.
                if var in self.derived_vars:
                    # Ex: {('PRECC', 'PRECL'): func, ('pr',): func1, ...}, is an OrderedDict.
                    possible_vars_and_funcs = self.derived_vars[var]

                    # Get the first valid variables and functions from possible vars.
                    # Ex: {('PRECC', 'PRECL'): func}
                    # These are checked to be in data_file.
                    vars_to_func_dict = self._get_first_valid_vars_climo(possible_vars_and_funcs, data_file)

                    # Get the variables as cdms2.TransientVariables.
                    # Ex: variables is [PRECC, PRECL], where both are cdms2.TransientVariables.
                    variables = self._get_original_vars_climo(vars_to_func_dict, data_file)

                    # Get the corresponding function.
                    # Ex: The func in {('PRECC', 'PRECL'): func}.
                    func = self._get_func(vars_to_func_dict)

                    # Call the function with the variables.
                    derived_var = func(*variables)

                # Or if the var is in the file, just get that.
                elif var in data_file.variables:
                    derived_var = data_file(var)(squeeze=1)

                # Otherwise, there's an error.
                else:
                    msg = 'Variable \'{}\' was not in the file {} nor was'.format(var, data_file.uri)
                    msg += ' it defined in the derived variables dictionary.'
                    raise RuntimeError(msg)

                variables.append(derived_var)

        return variables


    def _get_first_valid_vars_climo(self, vars_to_func_dict, data_file):
        """
        Given an OrderedDict of a list of variables to a function, ex: {('PRECC', 'PRECL'): func, (): func2},
        return the first valid {(vars): func} where the vars are in data_file.
        """
        vars_in_file = set(data_file.variables)
        possible_vars = vars_to_func_dict.keys()

        for list_of_vars in possible_vars:
            if vars_in_file.issuperset(list_of_vars):
                # All of the variables (list_of_vars) are in data_file.
                # Return the corresponding dict.
                return {list_of_vars: vars_to_func_dict[list_of_vars]}

        # In none of the derived variables are in the 
        msg = 'None of the variables in {} exist in the file {}'.format(possible_vars, data_file.uri)
        raise RuntimeError(msg)


    def _get_original_vars_climo(self, vars_to_func_dict, data_file):
        """
        Given a dictionary in the form {(vars): func}, get the vars
        from the data_file as cdms2.TransientVariables.

        These vars were checked to actually be in data_file.
        """
        # Since there's only one set of vars, we get the first
        # and only set of vars from the dictionary.
        vars_to_get = vars_to_func_dict.keys()[0]

        variables = [data_file(var)(squeeze=1) for var in vars_to_get]

        return variables


    def _get_func(self, vars_to_func_dict):
        """
        Get the function from the first and only entry in vars_to_func_dict,
        which is in the form {(vars): func}.
        """
        for k in vars_to_func_dict:
            return vars_to_func_dict[k]

 
    def _get_timeseries_var(self, data_path, start_yr, end_yr):
        """
        For a given season (self.season) and timeseries
        input data, get the variable (self.var).

        If self.extra_vars is also defined, get them as well.
        """
        pass

