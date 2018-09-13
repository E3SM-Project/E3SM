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
from . import general, climo


class Dataset():
    def __init__(self, parameters, ref=False, test=False, derived_vars={}, climo_fcn=None):
        self.parameters = parameters
        self.ref = ref
        self.test = test
        self.derived_vars = derived_vars
        self.climo_fcn = climo_fcn

        if self.ref is False and self.test is False:
            msg = 'Both ref and test cannot be False. One must be True.'
            raise RuntimeError(msg)
        elif self.ref is True and self.test is True:
            msg = 'Both ref and test cannot be True. Only one must be True.'
            raise RuntimeError(msg)

        if not self.derived_vars:
            # Use the default derived variables.
            self.derived_vars = acme_diags.derivations.acme.derived_variables

        if not self.climo_fcn:
            # Use the default climo function.
            self.climo_fcn = climo.climo

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
        self.extra_vars = extra_vars

        if not self.var:
            raise RuntimeError('Variable is invalid.')
        if not season:
            raise RuntimeError('Season is invalid.')

        # We need to make two decisions:
        # 1) Are the files being used reference or test data?
        #    - This is done with self.ref and self.test.
        # 2) Are the files being used climo or timeseries files?
        #   - This is done with the ref_timeseries_input and test_timeseries_input parameters.
        if self.ref and self.is_timeseries():
            # Get the reference variable from timeseries files.
            data_path = self.parameters.reference_data_path
            start_yr = self.parameters.ref_start_yr
            end_yr = self.parameters.ref_end_yr
            variables = self._get_timeseries_var(data_path, start_yr, end_yr, season)
            
        elif self.test and self.is_timeseries():
            # Get the test variable from timeseries files.
            data_path = self.parameters.test_data_path
            start_yr = self.parameters.test_start_yr
            end_yr = self.parameters.test_end_yr
            variables = self._get_timeseries_var(data_path, start_yr, end_yr, season)

        elif self.ref:
            # Get the reference variable from climo files.
            filename = general.get_ref_filename(self.parameters, season)
            variables = self._get_climo_var(filename)

        elif self.test:
            # Get the test variable from climo files.
            filename = general.get_test_filename(self.parameters, season)
            variables = self._get_climo_var(filename)

        else:
            msg = '''
            Error when determining what kind (ref or test) of variable to get and
            where to get it from (climo or timeseries files).
            '''
            raise RuntimeError(msg)

        # Needed so we can do:
        #   v1 = Dataset.get_variable('v1', season)
        # and also:
        #   v1, v2, v3 = Dataset.get_variable('v1', season, extra_vars=['v2', 'v3'])
        return variables[0] if len(variables) == 1 else variables


    def is_timeseries(self):
        """
        Return True if this dataset is for timeseries data.
        """
        if self.ref:
            return getattr(self.parameters, 'ref_timeseries_input', False)
        else:
            return getattr(self.parameters, 'test_timeseries_input', False)


    def is_climo(self):
        """
        Return True if this dataset is for climo data.
        """
        return not self.is_timeseries()


    def get_extra_variables_only(self, var, season, extra_vars):
        """
        For a given season, get only the extra variables.
        These can either be from the test data or reference data.
        """
        if not extra_vars:
            raise RuntimeError('Extra variables cannot be empty.')

        return self.get_variable(var, season, extra_vars)[1:]


    def get_attr_from_climo(self, attr, season):
        """
        For the given season, get the global attribute
        from the corresponding climo file.
        """
        if self.is_timeseries():
            raise RuntimeError('Cannot get a global attribute from timeseries files.')

        if self.ref:
            filename = general.get_ref_filename(self.parameters, season)
        else:
            filename = general.get_test_filename(self.parameters, season)
        
        with cdms2.open(filename) as f:
            return f.getglobal(attr)


    def get_start_end_time_slice(self):
        """
        Get the start and end years used when
        slicing the timeseries input data.
        """
        if not self.is_timeseries():
            msg = 'Input data is not timeseries, can\'t get start and end years.'
            raise RuntimeError(msg)

        # If neither *_start_time_slice or *_end_time_slice are defined,
        # then just used start_yr and end_yr.
        if self.ref:
            start_yr = self.parameters.ref_start_yr
            end_yr = self.parameters.ref_end_yr

            start_time_slice = getattr(self.parameters, 'ref_start_time_slice', start_yr)
            end_time_slice = getattr(self.parameters, 'ref_end_time_slice', end_yr)

        else:
            start_yr = self.parameters.test_start_yr
            end_yr = self.parameters.test_end_yr

            start_time_slice = getattr(self.parameters, 'test_start_time_slice', start_yr)
            end_time_slice = getattr(self.parameters, 'test_end_time_slice', end_yr)

        return start_time_slice, end_time_slice


    def _get_climo_var(self, filename):
        """
        For a given season and climo input data,
        get the variable (self.var).

        If self.extra_vars is also defined, get them as well.
        """
        vars_to_get = [self.var]
        vars_to_get.extend(self.extra_vars)
        return_variables = []

        with cdms2.open(filename) as data_file:
            for var in vars_to_get:
                # If it's a derived var, get that.
                if var in self.derived_vars:
                    # Ex: {('PRECC', 'PRECL'): func, ('pr',): func1, ...}, is an OrderedDict.
                    possible_vars_and_funcs = self.derived_vars[var]

                    # Get the first valid variables and functions from possible vars.
                    # Ex: {('PRECC', 'PRECL'): func}
                    # These are checked to be in data_file.
                    vars_to_func_dict = self._get_first_valid_vars_climo(possible_vars_and_funcs, data_file, var)

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
                    msg = 'Variable \'{}\' was not in the file {}, nor was'.format(var, data_file.uri)
                    msg += ' it defined in the derived variables dictionary.'
                    raise RuntimeError(msg)

                return_variables.append(derived_var)

        return return_variables


    def _get_first_valid_vars_climo(self, vars_to_func_dict, data_file, var):
        """
        Given an OrderedDict of a list of variables to a function
             ex: {('PRECC', 'PRECL'): func, ('var2',): func2},
        return the first valid {(vars): func} where the vars are in data_file.

        var is the actual variable the user requested.
        If none of the derived variables work, we try to just get this from the data_file.
        """
        vars_in_file = set(data_file.variables)
        possible_vars = vars_to_func_dict.keys()  # ex: [('pr',), ('PRECC', 'PRECL')]

        for list_of_vars in possible_vars:
            if vars_in_file.issuperset(list_of_vars):
                # All of the variables (list_of_vars) are in data_file.
                # Return the corresponding dict.
                return {list_of_vars: vars_to_func_dict[list_of_vars]}
        
        # None of the entries in the derived vars dictionary work,
        # so try to get the var directly.
        # Only try this if var actually exists in data_file.
        if var in data_file.variables:
            # The below will just cause var to get extracted from the data_file.
            return {(var,): lambda x: x}
        
        # Otherwise, there's no way to get the variable.
        msg = 'Neither does {} nor the variables in {}'.format(var, possible_vars)
        msg += ' exist in the file {}.'.format(data_file.uri)
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


    def _get_timeseries_var(self, data_path, start_yr, end_yr, season):
        """
        For a given season and timeseries input data,
        get the variable (self.var).

        If self.extra_vars is also defined, get them as well.
        """
        # Can't iterate through self.var and self.extra_vars as we do in _get_climo_var()
        # b/c the extra_vars must be taken from the same timeseries file as self.var.
        # So once we got a working vars_to_func_dict, we need to use this to get the extra_vars.

        return_variables = []

        # If it's a derived var, get that.
        if self.var in self.derived_vars:
            # Ex: {('PRECC', 'PRECL'): func, ('pr'): func1, ...}, is an OrderedDict.
            possible_vars_and_funcs = self.derived_vars[self.var]

            # Get the first valid variables and functions from possible vars.
            # Ex: {('PRECC', 'PRECL'): func}
            # These are checked, so there are valid timeseries files in data_path for these variables.
            vars_to_func_dict = self._get_first_valid_vars_timeseries(possible_vars_and_funcs, data_path, start_yr, end_yr)

            # Open the files of the variables and get the cdms2.TransientVariables.
            # Ex: [PRECC, PRECL], where both are TransientVariables.
            variables = self._get_original_vars_timeseries(vars_to_func_dict, data_path, start_yr, end_yr)

            # Get the corresponding function.
            # Ex: The func in {('PRECC', 'PRECL'): func}.
            func = self._get_func(vars_to_func_dict)

            # Call the function with the variables.
            derived_var = func(*variables)
            return_variables.append(derived_var)

            # Add any extra variables.
            # For a variable that is a derived variable, get all of the extra variables
            # from the 'first' original var.
            # Ex: We have {('PRECC', 'PRECL'): func} for PRECT.
            #     Any extra variables must come from PRECC_{start_yr}01_{end_yr}12.nc.
            first_orig_var = vars_to_func_dict.keys()[0][0]
            for extra_var in self.extra_vars:
                v = self._get_var_from_timeseries_file(first_orig_var, data_path, start_yr, end_yr, var_to_get=extra_var)
                return_variables.append(v)

        # Or if the timeseries file for the var exists, get that.
        elif self._does_timeseries_file_exist(self.var, data_path, start_yr, end_yr):
            # Find {var}_{start_yr}01_{end_yr}12.nc in data_path and get var from it.
            v = self._get_var_from_timeseries_file(self.var, data_path, start_yr, end_yr)
            return_variables.append(v)

            # Also get any extra vars.
            for extra_var in self.extra_vars:
                v = self._get_var_from_timeseries_file(self.var, data_path, start_yr, end_yr, var_to_get=extra_var)
                return_variables.append(v)

        # Otherwise, there's an error.
        else:
            msg = 'Variable \'{}\' doesn\'t have a file in the'.format(self.var)
            msg += ' directory {}, nor was'.format(data_path)
            msg += ' it defined in the derived variables dictionary.'
            raise RuntimeError(msg)

        # Call the climo on the variables.
        return [self.climo_fcn(v, season) for v in return_variables]


    def _get_first_valid_vars_timeseries(self, vars_to_func_dict, data_path, start_yr, end_yr):
        """
        Given an OrderedDict of a list of variables to a function
            ex: {('PRECC', 'PRECL'): func, ('var2',): func2},
        return the first valid {(vars): func} where the vars are variables from files in the form:
            {var}_{start_yr}01_{end_yr}12.nc
        located in data_path.

        If none of the derived variables work, we try to just get self.var in a file like:
            {self.var}_{start_yr}01_{end_yr}12.nc
        located in data_path.
        """
        possible_vars = vars_to_func_dict.keys()  # ex: [('pr',), ('PRECC', 'PRECL')]

        for list_of_vars in possible_vars:
            # Check that there are files in data_path that exist for all variables in list_of_vars.
            if all(self._does_timeseries_file_exist(var, data_path, start_yr, end_yr)
                    for var in list_of_vars):
                # All of the variables (list_of_vars) have files in data_path.
                # Return the corresponding dict.
                return {list_of_vars: vars_to_func_dict[list_of_vars]}
    
        # None of the entries in the derived vars dictionary are valid,
        # so try to get the var directly.
        # Only try this if there is a corresponding file for var in data_path.
        if self._does_timeseries_file_exist(self.var, data_path, start_yr, end_yr):
            # The below will just cause var to get extracted in {var}_{start_yr}01_{end_yr}12.nc.
            return {(self.var,): lambda x: x}

        # Otherwise, there's no way to get the variable.
        msg = 'Neither does {} nor the variables in {}'.format(self.var, possible_vars)
        msg += ' have valid files in {}.'.format(data_path)
        raise RuntimeError(msg)


    def _does_timeseries_file_exist(self, var, data_path, start_yr, end_yr):
        """
        Returns True if a file exists in data_path in the form:
            {var}_{start_yr}01_{end_yr}12.nc
        Or
            {self.parameters.ref_name}/{var}_{start_yr}01_{end_yr}12.nc
        """
        file_name = '{}_{}01_{}12.nc'.format(var, start_yr, end_yr)
        path = os.path.join(data_path, file_name)
        if os.path.exists(path):
            return True
        
        # Try looking for the file with the ref_name prepended to it.
        ref_name = getattr(self.parameters, 'ref_name', '')
        path = os.path.join(data_path, ref_name, file_name)
        # If it's not here, then the file doesn't exist.
        return os.path.exists(path)


    def _get_original_vars_timeseries(self, vars_to_func_dict, data_path, start_yr, end_yr):
        """
        Given a dictionary in the form {(vars): func}, get the vars
        from files in data_path as cdms2.TransientVariables.

        These vars were checked to actually be in data_path.
        start_yr and end_yr are used to get the filenames.
        """
        # Since there's only one set of vars, we get the first
        # and only set of vars from the dictionary.
        vars_to_get = vars_to_func_dict.keys()[0]

        variables = []
        for var in vars_to_get:
            v = self._get_var_from_timeseries_file(var, data_path, start_yr, end_yr)
            variables.append(v)
        
        return variables


    def _get_var_from_timeseries_file(self, var, data_path, start_yr, end_yr, var_to_get=''):
        """
        Get var from this file located in data_path:
            {var}_{start_yr}01_{end_yr}12.nc
        Or
            {self.parameters.ref_name}/{var}_{start_yr}01_{end_yr}12.nc
        If var_to_get is defined, get that from the file instead of var.
        """
        start_time_slice, end_time_slice = self.get_start_end_time_slice()

        start_time_slice = int(start_time_slice)
        end_time_slice = int(end_time_slice)
        
        file_name = '{}_{}01_{}12.nc'.format(var, start_yr, end_yr)
        path = os.path.join(data_path, file_name)
        if not os.path.exists(path):
            ref_name = getattr(self.parameters, 'ref_name', '')
            path = os.path.join(data_path, ref_name, file_name)

        var = var_to_get if var_to_get else var
        with cdms2.open(path) as f:
            return f(var, time=slice(start_time_slice, end_time_slice+1))(squeeze=1)
