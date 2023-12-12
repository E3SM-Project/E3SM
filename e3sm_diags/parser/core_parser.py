import argparse
import collections
import configparser
import copy
import hashlib
import itertools
import random
import re
import sys
import types
from io import StringIO
from typing import List, Optional

import yaml

from e3sm_diags.parameter.core_parameter import CoreParameter


class CoreParser:
    def __init__(self, parameter_cls=CoreParameter, *args, **kwargs):
        self.parser = argparse.ArgumentParser(  # type: ignore
            # "resolve" allows arguments to be overriden.
            conflict_handler="resolve",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            *args,
            **kwargs,
        )
        self.add_arguments()

        # The parser namespace, which is determined when parsing args.
        self._parser_namespace = None
        # The `e3sm_diags` parameter class associated with this instance.
        self._parameter_cls = parameter_cls

    @staticmethod
    def check_values_of_params(parameters: List[CoreParameter]):
        """
        Checks the parameter objects have all of the needed arguments with the
        correct values.

        This method loops over all parameter objects and calls their
        ``check_values()`` method.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of CoreParameter-based objects.
        """
        for p in parameters:
            p.check_values()

    def add_arguments(self):
        """Adds arguments to the parser object.

        A sub-class of `CoreParser` might extend this method with additional
        arguments unique to that sub-class or overwrite default arguments
        already defined in `CoreParser`.
        """
        self.parser.add_argument(
            "-p",
            "--parameters",
            type=str,
            dest="parameters",
            help="Path to the user-defined parameter file.",
            required=False,
        )
        self.parser.add_argument(
            "-d",
            "--diags",
            type=str,
            nargs="+",
            dest="other_parameters",
            default=[],
            help="Path to the other user-defined parameter file.",
            required=False,
        )
        self.parser.add_argument(
            "-n",
            "--num_workers",
            type=int,
            dest="num_workers",
            help="Number of workers, used when running with multiprocessing or in distributed mode.",
            required=False,
        )
        self.parser.add_argument(
            "--scheduler_addr",
            type=str,
            dest="scheduler_addr",
            help="Address of the scheduler in the form of IP_ADDRESS:PORT. Used when running in distributed mode.",
            required=False,
        )
        self.parser.add_argument(
            "-g",
            "--granulate",
            type=str,
            nargs="+",
            dest="granulate",
            help="A list of variables to granulate.",
            required=False,
        )
        self.parser.add_argument(
            "--selectors",
            type=str,
            nargs="+",
            dest="selectors",
            help="A list of variables to be used to select parameters from.",
            required=False,
        )

        self.parser.add_argument(
            "set_name",
            type=str,
            help="Name of the diags set to send " + "these arguments to.",
            nargs="?",
        )

        self.parser.add_argument(
            "-r",
            "--reference_data_set",
            type=str,
            dest="reference_data_set",
            help="List of observations or models that are used as a "
            + "reference against the test_data_set.",
            required=False,
        )

        self.parser.add_argument(
            "--reference_data_path",
            dest="reference_data_path",
            help="Path for the reference data.",
            required=False,
        )

        self.parser.add_argument(
            "--ref_timeseries_input",
            dest="ref_timeseries_input",
            help="The input reference data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--ref_start_yr",
            dest="ref_start_yr",
            help="Start year for the reference timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--ref_end_yr",
            dest="ref_end_yr",
            help="End year for the reference timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--ref_start_time_slice",
            dest="ref_start_time_slice",
            help="Starting time slice year for the reference timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--ref_end_time_slice",
            dest="ref_end_time_slice",
            help="Ending time slice year for the reference timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--ref_name",
            dest="ref_name",
            help="The string used to locate the reference file. "
            + "The reference file starts with the string.",
            required=False,
        )

        self.parser.add_argument(
            "--ref_file",
            dest="ref_file",
            help="Path to the reference file.",
            required=False,
        )

        self.parser.add_argument(
            "-t",
            "--test_data_set",
            type=str,
            dest="test_data_set",
            help="List of observations or models to test "
            + "against the reference_data_set.",
            required=False,
        )

        self.parser.add_argument(
            "--test_data_path",
            dest="test_data_path",
            help="Path for the test data.",
            required=False,
        )

        self.parser.add_argument(
            "--test_timeseries_input",
            dest="test_timeseries_input",
            help="The input test data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--test_start_yr",
            dest="test_start_yr",
            help="Start year for the test timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--test_end_yr",
            dest="test_end_yr",
            help="End year for the test timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--test_start_time_slice",
            dest="test_start_time_slice",
            help="Starting time slice year for the test timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--test_end_time_slice",
            dest="test_end_time_slice",
            help="Ending time slice year for the test timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--test_file",
            dest="test_file",
            help="Path to the test file.",
            required=False,
        )

        self.parser.add_argument(
            "--results_dir",
            dest="results_dir",
            help="Path of where to save the results.",
            required=False,
        )

        self.parser.add_argument(
            "--sets",
            nargs="+",
            dest="sets",
            help="Sets to use.",
            required=False,
        )

        self.parser.add_argument(
            "-D",
            "--dataset",
            dest="dataset",
            help="Dataset to use. Ex: 'ACME' or 'AMWG'.",
            required=False,
        )

        self.parser.add_argument(
            "--run_type",
            dest="run_type",
            help="What comparison to do. One of three options: "
            + "'model_vs_obs'/'obs_vs_model', 'model_vs_model', or 'obs_vs_obs'.",
            required=False,
        )

        self.parser.add_argument(
            "-v",
            "--variables",
            nargs="+",
            dest="variables",
            help="Variables to use.",
            required=False,
        )

        self.parser.add_argument(
            "--plevs",
            type=float,
            nargs="+",
            dest="plevs",
            help="Selected pressure level. [take list as input]",
            required=False,
        )

        self.parser.add_argument(
            "--plot_plevs",
            dest="plot_plevs",
            help="plot specified plevs",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--plot_log_plevs",
            dest="plot_log_plevs",
            help="plot plevs on log-scale",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "-s",
            "--seasons",
            nargs="+",
            dest="seasons",
            help="Seasons to use.",
            required=False,
        )

        self.parser.add_argument(
            "-r",
            "--regions",
            nargs="+",
            dest="regions",
            help="regions to use.",
            required=False,
        )

        self.parser.add_argument(
            "--regrid_tool",
            dest="regrid_tool",
            help="What regrid tool to use.",
            required=False,
        )

        self.parser.add_argument(
            "--regrid_method",
            dest="regrid_method",
            help="What regrid method for the regrid tool to use.",
            required=False,
        )

        self.parser.add_argument(
            "--case_id",
            dest="case_id",
            help="Defines a subdirectory to the metrics output, so multiple"
            + "cases can be compared.",
            required=False,
        )

        self.parser.add_argument(
            "--output_format",
            nargs="+",
            dest="output_format",
            help="What output format the plots should be saved in. "
            + "Possible values are: ['png', 'pdf', 'svg'].",
            required=False,
        )

        self.parser.add_argument(
            "--output_format_subplot",
            nargs="+",
            dest="output_format_subplot",
            help="What output format the individual subplots should be saved in (leave empty for no subplots)."
            + "Possible values are: ['png', 'pdf', 'svg'].",
            required=False,
        )

        self.parser.add_argument(
            "--canvas_size_w",
            type=int,
            dest="canvas_size_w",
            help="Size in pixels of the width for the output figure. " + "VCS only.",
            required=False,
        )

        self.parser.add_argument(
            "--canvas_size_h",
            type=int,
            dest="canvas_size_h",
            help="Size in pixels of the height for the output figure. " + "VCS only.",
            required=False,
        )

        self.parser.add_argument(
            "--figsize",
            type=float,
            nargs="+",
            dest="figsize",
            help="Width and height like so: [width, height]. " + "Matplotlib only.",
            required=False,
        )

        self.parser.add_argument(
            "--dpi",
            type=int,
            dest="dpi",
            help="DPI to use. " + "Matplotlib only.",
            required=False,
        )

        self.parser.add_argument(
            "--arrows",
            dest="arrows",
            help="Display arrows on the plot.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--logo",
            dest="logo",
            help="Display the logo. VCS only.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--contour_levels",
            type=float,
            nargs="+",
            dest="contour_levels",
            help="Levels for the test and reference plots.",
            required=False,
        )

        self.parser.add_argument(
            "--diff_levels",
            type=float,
            nargs="+",
            dest="diff_levels",
            help="Levels for the difference plot.",
            required=False,
        )

        self.parser.add_argument(
            "--reference_name",
            dest="reference_name",
            help="Name of the reference variable.",
            required=False,
        )

        self.parser.add_argument(
            "--test_name",
            dest="test_name",
            help="Name of the test variable.",
            required=False,
        )

        self.parser.add_argument(
            "--short_test_name",
            dest="short_test_name",
            help="User-defined test name.",
            required=False,
        )

        self.parser.add_argument(
            "--diff_name",
            dest="diff_name",
            help="Name of the difference variable.",
            required=False,
        )

        self.parser.add_argument(
            "--main_title",
            dest="main_title",
            help="The big title that appears on the top of the graph.",
            required=False,
        )

        self.parser.add_argument(
            "--reference_title",
            dest="reference_title",
            help="Title for the middle graph.",
            required=False,
        )

        self.parser.add_argument(
            "--test_title",
            dest="test_title",
            help="Title for the top graph.",
            required=False,
        )

        self.parser.add_argument(
            "--diff_title",
            dest="diff_title",
            help="Title for the bottom graph.",
            required=False,
        )

        self.parser.add_argument(
            "--reference_colormap",
            dest="reference_colormap",
            help="Colormap for the middle graph.",
            required=False,
        )

        self.parser.add_argument(
            "--test_colormap",
            dest="test_colormap",
            help="Colormap for the top graph.",
            required=False,
        )

        self.parser.add_argument(
            "--diff_colormap",
            dest="diff_colormap",
            help="Colormap for the bottom graph.",
            required=False,
        )

        self.parser.add_argument(
            "--reference_units",
            dest="reference_units",
            help="Units to use for the middle graph.",
            required=False,
        )

        self.parser.add_argument(
            "--test_units",
            dest="test_units",
            help="Units to use for the top graph.",
            required=False,
        )

        self.parser.add_argument(
            "--diff_units",
            dest="diff_units",
            help="Units to use for the bottom graph.",
            required=False,
        )

        self.parser.add_argument(
            "--backend",
            dest="backend",
            help="Graphical backend to use.",
            required=False,
        )

        self.parser.add_argument(
            "--multiprocessing",
            dest="multiprocessing",
            help="Run the diags using multiprocessing.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--save_netcdf",
            dest="save_netcdf",
            help="Save the NetCDF files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--no_viewer",
            dest="no_viewer",
            help="Don't generate the viewer.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--debug",
            dest="debug",
            help="Turns debugging on, allows code to prematurely break.",
            action="store_const",
            const=True,
            required=False,
        )

    def parse_args(
        self,
        args: Optional[List[str]] = None,
        namespace: Optional[argparse.Namespace] = None,
    ) -> argparse.Namespace:
        """Parses arguments passed from the command line.

        This method records the command used to run the parser, which is either
        `args` or `sys.argv`.

        Parameters
        ----------
        args : List[str], optional
            The list of arguments, by default [].
        namespace : Optional[argparse.Namespace], optional
            The namespace, by default None

        Returns
        -------
        Namespace
            The argparse.Namespace object with the parsed arguments.
        """
        # Remove arguments set by `ipykernel` via `sys.argv` because they
        # are not defined and recognized by `self.parser` using `add_argument()`.
        sys.argv = self._remove_ipykernel_args()
        self.cmd_used = sys.argv if args is None else args

        namespace, _ = self.parser.parse_known_args(args, namespace)

        return namespace

    def _remove_ipykernel_args(self) -> List[str]:
        """Removes `sys.argv` arguments set by `ipykernel`.

        `ipykernel` sets arguments during interactive console session, including
        when running scripts through Jupyter Notebook or using VSCode's
        interactive debugger. These arguments are not recognized by
        `self.parser`, which is an instance of `argparse.Argparser`. Since they
        are not recognized, the following error is raised:
        "ipykernel_launcher.py: error: unrecognized arguments: ...".

        Since `e3sm_diags` can be executed from the console via `e3sm_diags`
        command or a Python script, we need to drop these `ipykernel` arguments.
        Refer to [1]_ for more information.

        Returns
        -------
        List[str]
            A filtered list of args for `sys.argv`.

        References
        ----------
        .. [1] https://stackoverflow.com/a/67280166
        """
        regex_groups = [
            "(--ip=\\S*)",
            "(--stdin=\\S*)",
            "(--control=\\S*)",
            "(--hb=\\S*)",
            "(--Session.signature.scheme=\\S*)",
            "(--Session.key=\\S*)",
            "(--shell=\\S*)",
            "(--transport=\\S*)",
            "(--iopub=\\S*)",
            "(--ip=\\S*)",
            "(--f=\\S*\\/jupyter/\\S*.json)",
        ]
        pattern = ("|").join(regex_groups)
        regex = re.compile(pattern)

        new_sys_argv = [arg for arg in sys.argv if not regex.match(arg)]

        return new_sys_argv

    def view_args(self):
        """
        Returns the argparse.namespace
        """
        self._parse_arguments()
        return self._parser_namespace

    def get_parameters(
        self,
        cmdline_parameters=None,
        orig_parameters=None,
        other_parameters=[],
        default_vars=True,
        cmd_default_vars=True,
        *args,
        **kwargs,
    ):
        """
        Get the parameters based on the command line arguments and return a list of them.
        """
        if not cmdline_parameters:
            cmdline_parameters = self._get_cmdline_parameters(*args, **kwargs)
        if not orig_parameters:
            orig_parameters = self.get_orig_parameters(*args, **kwargs)
        if other_parameters == []:
            other_parameters = self.get_cfg_parameters(*args, **kwargs)

        # We don't want to add the selectors to each of the parameters.
        # Because if we do, it'll select all of the parameters at the end during the selection step.
        vars_to_ignore = self._get_selectors(
            cmdline_parameters, orig_parameters, other_parameters
        )
        self._combine_params(
            cmdline_parameters,
            orig_parameters,
            other_parameters,
            vars_to_ignore,
            default_vars,
            cmd_default_vars,
        )

        if other_parameters != []:
            final_parameters = other_parameters
        elif orig_parameters:
            final_parameters = [orig_parameters]
        elif cmdline_parameters:
            final_parameters = [cmdline_parameters]

        # User didn't give any command line options, so create a parameter from the
        # defaults of the command line argument or the Parameter class.
        elif cmd_default_vars:
            p = self._parameter_cls()
            for arg_name, arg_value in vars(self._parser_namespace).items():
                setattr(p, arg_name, arg_value)
            final_parameters = [p]
        elif default_vars:
            p = self._parameter_cls()
            final_parameters = [p]

        final_parameters = self._granulate(final_parameters)

        # Only select from the -p or the command line options.
        parameter = self.get_orig_parameters(*args, **kwargs)
        cmdline_parameter = self._get_cmdline_parameters(*args, **kwargs)

        # Command line parameters are added to parameter.
        # default_vars must be True, b/c the user excepts to select from them.
        self._combine_params(cmdline_parameter, parameter, default_vars=True)
        # Sometimes, one of these can be None, so get the one that's None.
        parameter = parameter if parameter else cmdline_parameter

        final_parameters = self.select(parameter, final_parameters)
        self._add_aliases(final_parameters)

        return final_parameters

    def get_orig_parameters(self, check_values=False, argparse_vals_only=True):
        """
        Returns the parameters created by -p. If -p wasn't used, returns None.
        """
        # FIXME: Should we deprecate this method for running `e3sm_diags`?
        self._parse_arguments()

        if not self._parser_namespace.parameters:  # type: ignore
            return None

        parameter = self._parameter_cls()

        # Remove all of the variables.
        parameter.__dict__.clear()

        # if self.__args_namespace.parameters is not None:
        parameter.load_parameter_from_py(self._parser_namespace.parameters)  # type: ignore

        if check_values:
            parameter.check_values()
        if argparse_vals_only:
            self._only_cmdline_args(parameter)

    def get_cfg_parameters(
        self,
        files_to_open: List[str] = [],
        check_values: bool = False,
        argparse_vals_only: bool = True,
    ) -> List[CoreParameter]:
        """Returns the parameters created using a diagnostic `cfg` file (-d).

        Diagnostic files should be `.cfg` format.

        Parameters
        ----------
        files_to_open : List[str], optional
            The diagnostic files to open, by default [].
            If ``files_to_open`` is set, then use the path specified instead of
            the file set by `-d`.
        check_values : bool, optional
            Check the values of the diagnostics file, by default False.
        argparse_vals_only : bool, optional
            Remove all parameters except those that are defined for the parser
            in ``self.add_arguments()``, by default True.

        Returns
        -------
        List[CoreParameter]
            A list of CoreParameter-based objects.

        Raises
        ------
        RuntimeError
            If the parameters input file is not `.json` or `.cfg` format.
        """

        parameters = []

        self._parse_arguments()

        if files_to_open == []:
            files_to_open = self._parser_namespace.other_parameters  # type: ignore

        if files_to_open is not None:
            for diags_file in files_to_open:
                if ".cfg" in diags_file:
                    params = self._get_cfg_parameters(
                        diags_file, check_values, argparse_vals_only
                    )
                else:
                    raise RuntimeError("The parameters input file must be a .cfg file")

                for p in params:
                    parameters.append(p)

        return parameters

    def _get_cfg_parameters(
        self, cfg_file, check_values=False, argparse_vals_only=True
    ):
        """
        Given a cfg file, return the parameters from it.
        """
        parameters = []

        cfg_file_obj = self._create_cfg_hash_titles(cfg_file)

        # Setting `strict=False` enables the parser to allow for any section
        # or duplicates while reading from a single source. This is required
        # because .cfg diagnostic files might contain duplicate sections with
        # slight tweaks based on the set.
        parser = configparser.ConfigParser(strict=False)
        parser.read_file(cfg_file_obj)

        for section in parser.sections():
            p = self._parameter_cls()

            # Remove all of the variables.
            p.__dict__.clear()

            for k, v in parser.items(section):
                v = yaml.safe_load(v)
                setattr(p, k, v)

            if check_values:
                p.check_values()
            if argparse_vals_only:
                self._only_cmdline_args(p)

            parameters.append(p)

        return parameters

    def _create_cfg_hash_titles(self, cfg_file):
        """
        Given a path to a cfg file, for any title '[#]', create a hash of it's contents
        and change the title to that. Then return the StringIO object.
        """
        lines = []
        with open(cfg_file) as f:
            lines = f.readlines()

        h_sha256 = hashlib.sha256()
        i = 0

        while i < len(lines):
            if lines[i] in ["[#]\n", "[#]"]:
                replace_idx = i
                str_list = []
                i += 1
                while i < len(lines) and not lines[i].startswith("["):
                    str_list.append(lines[i])
                    i += 1
                str_list.append(str(random.random()))  # Randomize the hash even more.
                h_sha256.update("".join(str_list).encode())
                lines[replace_idx] = "[{}]".format(h_sha256.hexdigest())
            else:
                i += 1
        return StringIO("\n".join(lines))

    def _get_cmdline_parameters(self, check_values=False, argparse_vals_only=True):
        """
        Use the other command line args besides -p and -d to create a single parameters object.
        """
        self._parse_arguments()

        if not self._were_cmdline_args_used():
            return None

        parameter = self._parameter_cls()

        # Remove all of the variables
        parameter.__dict__.clear()

        self._overwrite_parameters_with_cmdline_args(parameter)

        if check_values:
            parameter.check_values()
        if argparse_vals_only:
            self._only_cmdline_args(parameter)

        return parameter

    def _overwrite_parameters_with_cmdline_args(self, parameters):
        """
        Add the command line parameters used to the parameter object.
        """
        for arg_name, arg_value in vars(self._parser_namespace).items():
            if not self._is_arg_default_value(arg_name):
                setattr(parameters, arg_name, arg_value)

    def _were_cmdline_args_used(self):
        """
        Checks that other parameters, besides '-p' or '-d', were used.
        """
        for cmd in self.cmd_used:
            if cmd.startswith("-") and cmd not in [
                "-p",
                "--parameters",
                "-d",
                "--diags",
            ]:
                return True
        return False

    def _only_cmdline_args(self, parameter):
        """
        Remove all parameters except those that are usable by via command line
        arguments.
        """
        acceptable_args = vars(self.view_args())
        current_args = vars(parameter)

        params_to_del = [a for a in current_args if a not in acceptable_args]

        for param in params_to_del:
            delattr(parameter, param)

    def _parse_arguments(self):
        """
        Parse the command line arguments while checking for the user's arguments.
        """
        if self._parser_namespace is None:
            self._parser_namespace = self.parse_args()

    def select(self, main_parameters, parameters):
        """
        Given a list of parameters (parameters), only return those from this list
        whose 'selector' parameters are a subset of the 'selector' parameters of main_parameters.
        """

        def is_subset(param1, param2):
            """
            Check if param1 is a subset of param2.
            These are any Python objects.
            """
            if not isinstance(param1, list):
                param1 = [param1]
            if not isinstance(param2, list):
                param2 = [param2]

            return set(param1).issubset(set(param2))

        # Can't select from None.
        if not main_parameters:
            return parameters

        selectors = self._get_selectors(None, main_parameters, parameters)

        final_parameters = []

        for param in parameters:
            if all(
                is_subset(
                    getattr(param, select_parameter),
                    getattr(main_parameters, select_parameter),
                )
                for select_parameter in selectors
            ):
                final_parameters.append(param)

        return final_parameters

    def _get_selectors(
        self, cmdline_parameters=None, orig_parameters=None, other_parameters=None
    ):
        """
        Look through the cmdline_parameters, orig_parameters, and other_parameters
        in that order for the selectors used.
        If not defined in any of them, use the default one in the class.
        """
        if (
            hasattr(cmdline_parameters, "selectors")
            and cmdline_parameters.selectors is not None
        ):
            return cmdline_parameters.selectors
        elif (
            hasattr(orig_parameters, "selectors")
            and orig_parameters.selectors is not None
        ):
            return orig_parameters.selectors
        elif (
            hasattr(other_parameters, "selectors")
            and other_parameters.selectors is not None
        ):
            return other_parameters.selectors
        else:
            # If the parameter class has selectors, try to add that in.
            param = self._parameter_cls()
            if hasattr(param, "selectors"):
                return param.selectors
        # None of the passed in parameters have a selector and neither the main_parameter
        # nor the parameter class has selectors, so return an empty list.
        return []

    def _granulate(self, parameters):
        """
        Given a list of parameters objects, for each parameters with a `granulate` attribute,
        create multiple parameters objects for each result in the Cartesian product of `granulate`.
        """
        final_parameters = []
        for param in parameters:
            if not hasattr(param, "granulate") or (
                hasattr(param, "granulate") and not param.granulate
            ):
                final_parameters.append(param)
                continue

            # Remove any attrs that are modules from the param object.
            # These cause an error when copy.deepcopy(param) is used.
            attrs = vars(param).items()
            modules_in_param = []
            for var_name, var_value in attrs:
                if isinstance(var_value, types.ModuleType):
                    modules_in_param.append(var_name)
            for module in modules_in_param:
                delattr(param, module)

            # Granulate param.
            vars_to_granulate = param.granulate  # Ex: ['seasons', 'plevs']
            # Check that all of the vars_to_granulate are iterables.
            # Ex: {'season': ['ANN', 'DJF', 'MAM'], 'plevs': [850.0, 250.0]}
            vals_to_granulate = collections.OrderedDict()
            for v in vars_to_granulate:
                if not hasattr(param, v):
                    raise RuntimeError(
                        "Parameters object has no attribute '{}' to granulate.".format(
                            v
                        )
                    )
                param_v = getattr(param, v)
                if not isinstance(param_v, collections.abc.Iterable):
                    raise RuntimeError(
                        "Granulate option '{}' is not an iterable.".format(v)
                    )
                if param_v:  # Ignore [].
                    vals_to_granulate[v] = param_v

            # Ex: [('ANN', 850.0), ('ANN', 250.0), ('DJF', 850.0), ('DJF', 250.0), ...]
            granulate_values = list(itertools.product(*vals_to_granulate.values()))
            for g_vals in granulate_values:
                p = copy.deepcopy(param)
                for i, g_val in enumerate(g_vals):
                    key_at_index_i = list(vals_to_granulate.keys())[i]
                    # Make sure to insert a list with one element,
                    # which is why we have [g_val].
                    setattr(p, key_at_index_i, [g_val])
                final_parameters.append(p)

        return final_parameters

    def _combine_params(
        self,
        cmdline_parameters=None,
        orig_parameters=None,
        other_parameters=None,
        vars_to_ignore=[],
        default_vars=False,
        cmd_default_vars=False,
    ):
        """
        Combine cmdline_params (-* or --*), orig_parameters (-p), and other_parameters (-d),
        while ignoring any parameters listed in the 'selectors' parameter.
        Add any default arguments here as well.
        """
        if other_parameters:
            for parameters in other_parameters:
                self._add_default_values(parameters, default_vars, cmd_default_vars)

                # orig_parameters args take precedence over other_parameters.
                if orig_parameters:
                    for var in orig_parameters.__dict__:
                        if var not in vars_to_ignore:
                            parameters.__dict__[var] = orig_parameters.__dict__[var]

                # cmd_line args take the final precedence.
                if cmdline_parameters:
                    for var in cmdline_parameters.__dict__:
                        if var not in vars_to_ignore:
                            parameters.__dict__[var] = cmdline_parameters.__dict__[var]

        else:
            # Just combine cmdline_params with orig_params.
            if orig_parameters and cmdline_parameters:
                self._add_default_values(
                    orig_parameters, default_vars, cmd_default_vars
                )

                for var in cmdline_parameters.__dict__:
                    # if var not in vars_to_ignore and self._was_command_used(var):
                    if var not in vars_to_ignore:
                        # Only add it if it was not in param and was passed from cmd line.
                        orig_parameters.__dict__[var] = cmdline_parameters.__dict__[var]

            elif orig_parameters:
                self._add_default_values(
                    orig_parameters, default_vars, cmd_default_vars
                )
            elif cmdline_parameters:
                self._add_default_values(
                    cmdline_parameters, default_vars, cmd_default_vars
                )

    def _add_default_values(
        self, parameter, default_vars=False, cmd_default_vars=False
    ):
        """
        Add the default values to the parameter.
        These can come from the default values defined in the Parameter class,
        or the `default` option defined in ArgumentParser.add_argument().
        """
        # Add the command line default parameters first.
        if cmd_default_vars:
            for arg_name, arg_value in vars(self._parser_namespace).items():
                if arg_name in parameter.__dict__ or not self._is_arg_default_value(
                    arg_name
                ):
                    continue
                # Only add the default values, that aren't already in parameter.
                setattr(parameter, arg_name, arg_value)

        # Then add the defaults defined in the Parameter class.
        if default_vars:
            for arg_name, arg_value in vars(self._parameter_cls()).items():
                if arg_name in parameter.__dict__:
                    continue
                setattr(parameter, arg_name, arg_value)

    def _add_aliases(self, parameters):
        """
        For each of the parameters, add all of
        the defined aliases as other attributes.
        """
        for param in parameters:
            # We need to set this info as a variable.
            # Can't do:
            #    for param_name in vars(param)
            # because we're modifying param as we iterate.
            # We also need to make a copy because dicts are referenced.
            param_names = copy.copy(vars(param))

            for param_name in param_names:
                param_value = getattr(param, param_name)
                aliases = self._get_alias(param_name)
                # Add all of the aliases for param_name to the param object.
                for alias in aliases:
                    setattr(param, alias, param_value)

    def _get_alias(self, param):
        """
        For a single parameter, get the aliases of it.
        """
        # Parameters can start with either '-' or '--'.
        param = "--{}".format(param)
        if param not in self.parser._option_string_actions:
            param = "-{}".format(param)
        if param not in self.parser._option_string_actions:
            return []

        # Ex: If param is 'parameters', then we get ['-p', '--parameters'].
        aliases = self.parser._option_string_actions[param].option_strings

        return [a.replace("-", "") for a in aliases]

    def _is_arg_default_value(self, arg):
        """
        Look at the command used for this parser (ex: test.py -s something --s1 something1)
        and if arg wasn't used, then it's a default value.
        """
        # Each cmdline_arg is either '-*' or '--*'.
        for cmdline_arg in self.parser._option_string_actions:
            if arg == self.parser._option_string_actions[
                cmdline_arg
            ].dest and self._was_command_used(cmdline_arg):
                return False
        return True

    def _was_command_used(self, cmdline_arg):
        """
        Returns True if the cmdline_arg was used to
        run the script that has this parser.
        """
        # self.cmd_used is like: ['something.py', '-p', 'test.py', '--s1', 'something']
        for cmd in self.cmd_used:
            # Sometimes, a command is run with '=': 'driver.py --something=this'
            for c in cmd.split("="):
                if cmdline_arg == c:
                    return True
        return False
