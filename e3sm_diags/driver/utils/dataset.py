"""
Get a variable from input data (either reference or test data).
This data can either be climatology files or timeseries files.
Derived variables are also supported.
"""
import collections
import fnmatch
import glob
import os
import re

import cdms2

import e3sm_diags.derivations.acme
from e3sm_diags.driver import utils

from . import climo


class Dataset:
    def __init__(
        self,
        parameters,
        ref=False,
        test=False,
        derived_vars={},
        climo_fcn=None,
    ):
        self.parameters = parameters
        self.ref = ref
        self.test = test
        self.derived_vars = derived_vars
        self.climo_fcn = climo_fcn

        if self.ref is False and self.test is False:
            msg = "Both ref and test cannot be False. One must be True."
            raise RuntimeError(msg)
        elif self.ref is True and self.test is True:
            msg = "Both ref and test cannot be True. Only one must be True."
            raise RuntimeError(msg)

        if not self.derived_vars:
            # Use the default derived variables.
            self.derived_vars = e3sm_diags.derivations.acme.derived_variables

        if not self.climo_fcn:
            # Use the default climo function.
            self.climo_fcn = climo.climo

        if hasattr(self.parameters, "derived_variables"):
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

    def get_timeseries_variable(
        self, var, extra_vars=[], single_point=False, *args, **kwargs
    ):
        """
        Get the variable and any extra variables, only if they are timeseries files.
        These variables can either be from the test data or reference data.
        """
        self.var = var
        self.extra_vars = extra_vars

        if not self.is_timeseries():
            msg = "You can only use this function with timeseries data."
            raise RuntimeError(msg)

        if self.ref:
            # Get the reference variable from timeseries files.
            data_path = self.parameters.reference_data_path
            variables = self._get_timeseries_var(data_path, *args, **kwargs)

        elif self.test:
            # Get the test variable from timeseries files.
            data_path = self.parameters.test_data_path
            variables = self._get_timeseries_var(data_path, *args, **kwargs)

        else:
            msg = "Error when determining what kind (ref or test)of variable to get."
            raise RuntimeError(msg)

        # Needed so we can do:
        #   v1 = Dataset.get_variable('v1', season)
        # and also:
        #   v1, v2, v3 = Dataset.get_variable('v1', season, extra_vars=['v2', 'v3'])

        # Need to double check sub_monthly flag when applying to sub_monthly time series later
        sub_monthly = False

        if single_point:
            sub_monthly = True

        for variable in variables:
            if variable.getTime() and not sub_monthly:
                variable = utils.general.adjust_time_from_time_bounds(variable)
        return variables[0] if len(variables) == 1 else variables

    def get_climo_variable(self, var, season, extra_vars=[], *args, **kwargs):
        """
        For a given season, get the variable and any extra variables and run
        the climatology on them.
        These variables can either be from the test data or reference data.
        """
        self.var = var
        self.extra_vars = extra_vars

        if not self.var:
            raise RuntimeError("Variable is invalid.")
        if not season:
            raise RuntimeError("Season is invalid.")

        # We need to make two decisions:
        # 1) Are the files being used reference or test data?
        #    - This is done with self.ref and self.test.
        # 2) Are the files being used climo or timeseries files?
        #   - This is done with the ref_timeseries_input and test_timeseries_input parameters.
        if self.ref and self.is_timeseries():
            # Get the reference variable from timeseries files.
            data_path = self.parameters.reference_data_path
            timeseries_vars = self._get_timeseries_var(data_path, *args, **kwargs)
            # Run climo on the variables.
            variables = [self.climo_fcn(v, season) for v in timeseries_vars]

        elif self.test and self.is_timeseries():
            # Get the test variable from timeseries files.
            data_path = self.parameters.test_data_path
            timeseries_vars = self._get_timeseries_var(data_path, *args, **kwargs)
            # Run climo on the variables.
            variables = [self.climo_fcn(v, season) for v in timeseries_vars]

        elif self.ref:
            # Get the reference variable from climo files.
            filename = self.get_ref_filename_climo(season)
            variables = self._get_climo_var(filename, *args, **kwargs)

        elif self.test:
            # Get the test variable from climo files.
            filename = self.get_test_filename_climo(season)
            variables = self._get_climo_var(filename, *args, **kwargs)

        else:
            msg = "Error when determining what kind (ref or test) "
            msg += "of variable to get and where to get it from "
            msg += "(climo or timeseries files)."
            raise RuntimeError(msg)

        # Needed so we can do:
        #   v1 = Dataset.get_variable('v1', season)
        # and also:
        #   v1, v2, v3 = Dataset.get_variable('v1', season, extra_vars=['v2', 'v3'])
        return variables[0] if len(variables) == 1 else variables

    def get_static_variable(self, static_var, primary_var):
        if self.ref:
            # Get the reference variable from timeseries files.
            data_path = self.parameters.reference_data_path
        elif self.test:
            # Get the test variable from timeseries files.
            data_path = self.parameters.test_data_path
        file_path = self._get_timeseries_file_path(primary_var, data_path)
        fin = cdms2.open(file_path)
        result = fin(static_var)
        fin.close()
        return result

    def is_timeseries(self):
        """
        Return True if this dataset is for timeseries data.
        """
        if self.ref:
            return getattr(self.parameters, "ref_timeseries_input", False)
        else:
            return getattr(self.parameters, "test_timeseries_input", False)

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
            raise RuntimeError("Extra variables cannot be empty.")

        if self.is_climo():
            return self.get_climo_variable(
                var, season, extra_vars, extra_vars_only=True
            )
        else:
            return self.get_timeseries_variable(var, extra_vars, extra_vars_only=True)

        return self.get_climo_variable(var, season, extra_vars, extra_vars_only=True)

    def get_attr_from_climo(self, attr, season):
        """
        For the given season, get the global attribute
        from the corresponding climo file.
        """
        if self.is_timeseries():
            raise RuntimeError("Cannot get a global attribute from timeseries files.")

        if self.ref:
            filename = self.get_ref_filename_climo(season)
        else:
            filename = self.get_test_filename_climo(season)

        with cdms2.open(filename) as f:
            return f.getglobal(attr)

    def get_start_and_end_years(self):
        """
        Get the user-defined start and end years.
        """
        sub_monthly = False
        if self.parameters.sets[0] in ["area_mean_time_series"]:
            start_yr = getattr(self.parameters, "start_yr")
            end_yr = getattr(self.parameters, "end_yr")
        else:
            if self.ref:
                start_yr = getattr(self.parameters, "ref_start_yr")
                end_yr = getattr(self.parameters, "ref_end_yr")

            else:
                start_yr = getattr(self.parameters, "test_start_yr")
                end_yr = getattr(self.parameters, "test_end_yr")

        if self.parameters.sets[0] in ["diurnal_cycle", "arm_diags"]:
            sub_monthly = True

        return start_yr, end_yr, sub_monthly

    def get_test_filename_climo(self, season):
        """
        Return the path to the test file name based on
        the season and other parameters.
        For climo files only.
        """
        path = self.parameters.test_data_path
        data_name = self.parameters.test_name

        if hasattr(self.parameters, "test_file"):
            fnm = os.path.join(path, self.parameters.test_file)
            if not os.path.exists(fnm):
                raise IOError("File not found: {}".format(fnm))
            return fnm

        return self._get_climo_filename(path, data_name, season)

    def get_ref_filename_climo(self, season):
        """
        Return the path to the reference file name based on
        the season and other parameters.
        For climo files only.
        """
        path = self.parameters.reference_data_path
        data_name = self.parameters.ref_name

        if (
            hasattr(self.parameters, "ref_file")
            and getattr(self.parameters, "ref_file") != ""
        ):
            fnm = os.path.join(path, self.parameters.ref_file)
            if not os.path.exists(fnm):
                raise IOError("File not found: {}".format(fnm))
            return fnm

        return self._get_climo_filename(path, data_name, season)

    def _get_climo_filename(self, path, data_name, season):
        """
        For climo files, return the path of the file based on the parameters.
        If the file isn't found, try looking for it in path/data_name/ dir as well.
        """
        fnm = self._find_climo_file(path, data_name, season)
        if not os.path.exists(fnm):
            # Try looking for the file nested in a folder, based on the test_name.
            pth = os.path.join(path, data_name)
            if os.path.exists(pth):
                fnm = self._find_climo_file(pth, data_name, season)

        if not os.path.exists(fnm):
            raise IOError(
                "No file found for {} and {} in {}".format(data_name, season, path)
            )

        return fnm

    def _find_climo_file(self, path_name, data_name, season):
        """
        Locate climatology file name based on data_name and season.
        """
        dir_files = sorted(os.listdir(path_name))
        for filename in dir_files:
            if filename.startswith(data_name + "_" + season):
                return os.path.join(path_name, filename)
        # The below is only ran on model data, because a shorter name is passed into this software. Won't work when use month name such as '01' as season.
        for filename in dir_files:
            if season in [
                "ANN",
                "DJF",
                "MAM",
                "JJA",
                "SON",
                "01",
                "02",
                "03",
                "04",
                "05",
                "06",
                "07",
                "08",
                "09",
                "10",
                "11",
                "12",
            ]:
                if filename.startswith(data_name) and season in filename:
                    return os.path.join(path_name, filename)
        # No file found.
        return ""

    def _get_climo_var(self, filename, extra_vars_only=False):
        """
        For a given season and climo input data,
        get the variable (self.var).

        If self.extra_vars is also defined, get them as well.
        """
        vars_to_get = []
        if not extra_vars_only:
            vars_to_get.append(self.var)
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
                    vars_to_func_dict = self._get_first_valid_vars_climo(
                        possible_vars_and_funcs, data_file, var
                    )

                    # Get the variables as cdms2.TransientVariables.
                    # Ex: variables is [PRECC, PRECL], where both are cdms2.TransientVariables.
                    variables = self._get_original_vars_climo(
                        vars_to_func_dict, data_file
                    )

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
                    msg = "Variable '{}' was not in the file {}, nor was".format(
                        var, data_file.uri
                    )
                    msg += " it defined in the derived variables dictionary."
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
        possible_vars = list(
            vars_to_func_dict.keys()
        )  # ex: [('pr',), ('PRECC', 'PRECL')]

        # Add support for wild card `?` in variable strings: ex ('bc_a?DDF', 'bc_c?DDF')
        for list_of_vars in possible_vars:
            matched_var_list = list(list_of_vars).copy()
            for var_list in list_of_vars:
                if "?" in var_list:
                    matched_var_list += fnmatch.filter(list(vars_in_file), var_list)
                    matched_var_list.remove(var_list)

            if vars_in_file.issuperset(tuple(matched_var_list)):
                # All of the variables (list_of_vars) are in data_file.
                # Return the corresponding dict.
                return {tuple(matched_var_list): vars_to_func_dict[list_of_vars]}

        # None of the entries in the derived vars dictionary work,
        # so try to get the var directly.
        # Only try this if var actually exists in data_file.
        if var in data_file.variables:
            # The below will just cause var to get extracted from the data_file.
            return {(var,): lambda x: x}

        # Otherwise, there's no way to get the variable.
        msg = "Neither does {} nor the variables in {}".format(var, possible_vars)
        msg += " exist in the file {}.".format(data_file.uri)
        raise RuntimeError(msg)

    def _get_original_vars_climo(self, vars_to_func_dict, data_file):
        """
        Given a dictionary in the form {(vars): func}, get the vars
        from the data_file as cdms2.TransientVariables.

        These vars were checked to actually be in data_file.
        """
        # Since there's only one set of vars, we get the first
        # and only set of vars from the dictionary.
        vars_to_get = list(vars_to_func_dict.keys())[0]

        variables = [data_file(var)(squeeze=1) for var in vars_to_get]

        return variables

    def _get_func(self, vars_to_func_dict):
        """
        Get the function from the first and only entry in vars_to_func_dict,
        which is in the form {(vars): func}.
        """
        for k in vars_to_func_dict:
            return vars_to_func_dict[k]

    def _get_timeseries_var(self, data_path, extra_vars_only=False):
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
            vars_to_func_dict = self._get_first_valid_vars_timeseries(
                possible_vars_and_funcs, data_path
            )

            # We do want the self.var.
            if not extra_vars_only:
                # Open the files of the variables and get the cdms2.TransientVariables.
                # Ex: [PRECC, PRECL], where both are TransientVariables.
                variables = self._get_original_vars_timeseries(
                    vars_to_func_dict, data_path
                )

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
            first_orig_var = list(vars_to_func_dict.keys())[0][0]
            for extra_var in self.extra_vars:
                v = self._get_var_from_timeseries_file(
                    first_orig_var, data_path, var_to_get=extra_var
                )
                return_variables.append(v)

        # Or if the timeseries file for the var exists, get that.
        elif self._get_timeseries_file_path(self.var, data_path):
            # We do want the self.var.
            if not extra_vars_only:
                # Find {var}_{start_yr}01_{end_yr}12.nc in data_path and get var from it.
                v = self._get_var_from_timeseries_file(self.var, data_path)
                return_variables.append(v)

            # Also get any extra vars.
            for extra_var in self.extra_vars:
                v = self._get_var_from_timeseries_file(
                    self.var, data_path, var_to_get=extra_var
                )
                return_variables.append(v)

        # Otherwise, there's an error.
        else:
            msg = "Variable '{}' doesn't have a file in the".format(self.var)
            msg += " directory {}, nor was".format(data_path)
            msg += " it defined in the derived variables dictionary."
            raise RuntimeError(msg)

        return return_variables

    def _get_first_valid_vars_timeseries(self, vars_to_func_dict, data_path):
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
        possible_vars = list(
            vars_to_func_dict.keys()
        )  # ex: [('pr',), ('PRECC', 'PRECL')]

        for list_of_vars in possible_vars:
            # Check that there are files in data_path that exist for all variables in list_of_vars.
            if all(
                self._get_timeseries_file_path(var, data_path) for var in list_of_vars
            ):
                # All of the variables (list_of_vars) have files in data_path.
                # Return the corresponding dict.
                return {list_of_vars: vars_to_func_dict[list_of_vars]}

        # None of the entries in the derived vars dictionary are valid,
        # so try to get the var directly.
        # Only try this if there is a corresponding file for var in data_path.
        if self._get_timeseries_file_path(self.var, data_path):
            # The below will just cause var to get extracted in {var}_{start_yr}01_{end_yr}12.nc.
            return {(self.var,): lambda x: x}

        # Otherwise, there's no way to get the variable.
        msg = "Neither does {} nor the variables in {}".format(self.var, possible_vars)
        msg += " have valid files in {}.".format(data_path)
        raise RuntimeError(msg)

    def _get_timeseries_file_path(self, var, data_path):
        """
        Returns the file path if a file exists in data_path in the form:
            {var}_{start_yr}01_{end_yr}12.nc
        Or
            {self.parameters.ref_name}/{var}_{start_yr}01_{end_yr}12.nc
        This is equivalent to returning True if the file exists.

        If there are multiple files that exist for a variable
        (with different start_yr or end_yr), return ''.
        This is equivalent to returning False.
        """
        # Get all of the nc file paths in data_path.

        # path = os.path.join(data_path, '*.nc')
        path = os.path.join(data_path, "*.*")
        files = sorted(glob.glob(path))

        # Both .nc and .xml files are supported
        file_fmt = ""
        if len(files) > 0:
            file_fmt = files[0].split(".")[-1]

        # Everything between '{var}_' and '.nc' in a
        # time-series file is always 13 characters.
        if self.parameters.sets[0] in ["arm_diags"]:
            site = getattr(self.parameters, "regions", "")
            re_str = var + "_" + site[0] + r"_.{13}." + file_fmt
        else:
            re_str = var + r"_.{13}." + file_fmt
        re_str = os.path.join(data_path, re_str)
        matches = [f for f in files if re.search(re_str, f)]

        if len(matches) == 1:
            return matches[0]
        elif len(matches) >= 2:
            msg = "For the variable {} you have two timeseries files in the ".format(
                var
            )
            msg += "directory: {} This currently isn't supported.".format(data_path)
            raise RuntimeError(msg)

        # If nothing was found, try looking for the file with
        # the ref_name prepended to it.
        ref_name = getattr(self.parameters, "ref_name", "")
        # path = os.path.join(data_path, ref_name, '*.nc')
        path = os.path.join(data_path, ref_name, "*.*")
        files = sorted(glob.glob(path))
        # Both .nc and .xml files are supported
        file_fmt = ""
        if len(files) > 0:
            file_fmt = files[0].split(".")[-1]

        # Everything between '{var}_' and '.nc' in a
        # time-series file is always 13 characters.
        re_str = var + r"_.{13}." + file_fmt
        re_str = os.path.join(data_path, ref_name, re_str)
        matches = [f for f in files if re.search(re_str, f)]
        # Again, there should only be one file per var in this new location.
        if len(matches) == 1:
            return matches[0]
        elif len(matches) >= 2:
            msg = "For the variable {} you have two timeseries files in the ".format(
                var
            )
            msg += "directory: {} This currently isn't supported.".format(data_path)
            raise RuntimeError(msg)
        else:
            return ""

    def _get_original_vars_timeseries(self, vars_to_func_dict, data_path):
        """
        Given a dictionary in the form {(vars): func}, get the vars
        from files in data_path as cdms2.TransientVariables.

        These vars were checked to actually be in
        data_path in _get_first_valid_vars_timeseries().
        """
        # Since there's only one set of vars, we get the first
        # and only set of vars from the dictionary.
        vars_to_get = list(vars_to_func_dict.keys())[0]

        variables = []
        for var in vars_to_get:
            v = self._get_var_from_timeseries_file(var, data_path)
            variables.append(v)

        return variables

    def _get_var_from_timeseries_file(self, var, data_path, var_to_get=""):
        """
        Get the actual var from the timeseries file for var.
        If var_to_get is defined, get that from the file instead of var.

        This function is only called after it's checked that a file
        for this var exists in data_path.
        The checking is done in _get_first_valid_vars_timeseries().
        """
        (
            start_year,
            end_year,
            sub_monthly,
        ) = self.get_start_and_end_years()
        if sub_monthly:
            start_time = "{}-01-01".format(start_year)
            end_time = "{}-01-01".format(str(int(end_year) + 1))
            slice_flag = "co"
        else:
            start_time = "{}-01-15".format(start_year)
            end_time = "{}-12-15".format(end_year)
            slice_flag = "ccb"

        fnm = self._get_timeseries_file_path(var, data_path)

        var = var_to_get if var_to_get else var

        # get available start and end years from file name: {var}_{start_yr}01_{end_yr}12.nc
        start_year = int(start_year)
        end_year = int(end_year)
        var_start_year = int(fnm.split("/")[-1].split("_")[-2][:4])
        var_end_year = int(fnm.split("/")[-1].split("_")[-1][:4])

        if start_year < var_start_year:
            msg = "Invalid year range specified for test/reference time series data: start_year={}<{}=var_start_yr".format(
                start_year, var_start_year
            )
            raise RuntimeError(msg)
        elif end_year > var_end_year:
            msg = "Invalid year range specified for test/reference time series data: end_year={}>{}=var_end_yr".format(
                end_year, var_end_year
            )
            raise RuntimeError(msg)
        else:
            # with cdms2.open(fnm) as f:
            #    var_time = f(var, time=(start_time, end_time, 'ccb'))(squeeze=1)
            #    return var_time
            # For xml files using above with statement won't work because the Dataset object returned doesn't have attribute __enter__ for content management.
            fin = cdms2.open(fnm)
            var_time = fin(var, time=(start_time, end_time, slice_flag))(squeeze=1)
            fin.close()
            return var_time
