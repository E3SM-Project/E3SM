import copy
from itertools import chain
from typing import List, Union

import e3sm_diags  # noqa: F401
from e3sm_diags.e3sm_diags_driver import get_default_diags_path, main
from e3sm_diags.logger import custom_logger, move_log_to_prov_dir
from e3sm_diags.parameter import SET_TO_PARAMETERS
from e3sm_diags.parameter.core_parameter import DEFAULT_SETS, CoreParameter
from e3sm_diags.parser.core_parser import CoreParser

logger = custom_logger(__name__)


class Run:
    """
    Used to run diagnostics.
    A class is needed because we often need to store some
    state regarding what we need to run, like the sets selected.
    """

    def __init__(self):
        self.parser = CoreParser()

        # The list of sets to run based on diagnostic parameters.
        self.sets_to_run = []

        # The path to the user-specified `.cfg` file using `-d/--diags` or
        # the default diagnostics `.cfg` file.
        self.cfg_path = None

    @property
    def is_cfg_file_arg_set(self):
        """A property to check if `-d/--diags` was set to a `.cfg` filepath.

        Returns
        -------
        bool
            True if list contains more than one path, else False.
        """
        args = self.parser.view_args()
        self.cfg_path = args.other_parameters

        is_set = len(self.cfg_path) > 0

        return is_set

    def run_diags(
        self, parameters: List[CoreParameter], use_cfg: bool = True
    ) -> Union[List[CoreParameter], None]:
        """Run a set of diagnostics with a list of parameters.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of parameters defined through the Python API.
        use_cfg : bool, optional
            Also run diagnostics using a `.cfg` file, by default True.

              * If True, run all sets using the list of parameters passed in
                this function and parameters defined in a .cfg file (if
                defined), or use the .cfg file(s) for default diagnostics. This
                is the default option.
              * If False, only the parameters passed via ``parameters`` will be
                run. The sets to run are based on the sets defined by the
                parameters. This makes it easy to debug a few sets instead of
                all of the debug sets too.

        Returns
        -------
        Union[List[CoreParameter], None]
            A list of parameter objects with their results (if successful).

        Raises
        ------
        RuntimeError
            If a diagnostic run using a parameter fails for any reason.
        """
        params = self.get_run_parameters(parameters, use_cfg)
        params_results = None

        if params is None or len(params) == 0:
            raise RuntimeError(
                "No parameters we able to be extracted. Please "
                "check the parameters you defined."
            )

        try:
            params_results = main(params)
        except Exception:
            logger.exception("Error traceback:", exc_info=True)

        # param_results might be None because the run(s) failed, so move
        # the log using the `params[0].results_dir` instead.
        move_log_to_prov_dir(params[0].results_dir)

        return params_results

    def get_run_parameters(
        self, parameters: List[CoreParameter], use_cfg: bool = True
    ) -> List[CoreParameter]:
        """Get the run parameters.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of parameters defined through the Python API.
        use_cfg : bool, optional
            Use parameters defined in a .cfg file too, by default True.

        Returns
        -------
        List[CoreParameter]
            A list of run parameters.
        """
        self._validate_parameters(parameters)

        # FIXME: This line produces some unintended side-effects. For example,
        # let's say we have two objects: 1. CoreParameter, 2. ZonalMean2DParameter.
        # If object 1 has `plevs=[200]`, this value will get copied to object 2.
        # Object 2 has a check to make sure plevs has more than 1 value. This
        # breaks the diagnostic run as a result. The workaround is to loop
        # over `run_diags()` function and run one parameter at a time.
        self._add_parent_attrs_to_children(parameters)

        if use_cfg:
            run_params = self._get_parameters_with_api_and_cfg(parameters)
        else:
            run_params = self._get_parameters_with_api_only(parameters)

        self.parser.check_values_of_params(run_params)

        return run_params

    def _validate_parameters(self, parameters: List[CoreParameter]):
        if parameters is None or not isinstance(parameters, list):
            raise RuntimeError("You must pass in a list of parameter objects.")

        param_types_list = [
            p.__class__ for p in parameters if p.__class__ != CoreParameter
        ]
        param_types_set = set(param_types_list)

        if len(param_types_set) != len(param_types_list):
            raise RuntimeError(
                "You passed in two or more non-CoreParameter objects of the same type."
            )

    def _get_parameters_with_api_and_cfg(
        self, parameters: List[CoreParameter]
    ) -> List[CoreParameter]:
        """Get the run parameters using the Python API and .cfg file.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of parameter objects defined by the Python API.

        Returns
        -------
        List[CoreParameter]
            A list of run parameters, including ones defined in a .cfg file.
            Any non-CoreParameter objects will be replaced by a sub-class
            based on the set (``SETS_TO_PARAMETERS``).
        """
        run_params = []

        if len(self.sets_to_run) == 0:
            self.sets_to_run = DEFAULT_SETS

        for set_name in self.sets_to_run:
            if self.is_cfg_file_arg_set:
                cfg_params = self._get_custom_params_from_cfg_file()
            else:
                run_type = parameters[0].run_type
                cfg_params = self._get_default_params_from_cfg_file(run_type)

            param = self._get_instance_of_param_class(
                SET_TO_PARAMETERS[set_name], parameters
            )

            # Since each parameter will have lots of default values, we want to
            # remove them. Otherwise when calling get_parameters(), these
            # default values will take precedence over values defined in
            # other_params.
            self._remove_attrs_with_default_values(param)
            param.sets = [set_name]

            params = self.parser.get_parameters(
                orig_parameters=param,
                other_parameters=cfg_params,
                cmd_default_vars=False,
                argparse_vals_only=False,
            )

            # Makes sure that any parameters that are selectors will be in param.
            self._add_attrs_with_default_values(param)

            # The select() call in get_parameters() was made for the original
            # command-line way of using CDP. We just call it manually with the
            # parameter object param.
            params = self.parser.select(param, params)

            run_params.extend(params)

        return run_params

    def _get_custom_params_from_cfg_file(self) -> List[CoreParameter]:
        """Get parameters using the cfg file set by `-d`/`--diags`.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects.
        """
        params = self.parser.get_cfg_parameters(argparse_vals_only=False)

        params_final = self._convert_params_to_subclass(params)

        return params_final

    def _get_default_params_from_cfg_file(self, run_type: str) -> List[CoreParameter]:
        """Get parameters using the default diagnostic .cfg file(s).

        Parameters
        ----------
        run_type : str
            The run type used to check for which .cfg file(s) to reference.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects.
        """
        # Get the paths to the default diags .cfg file(s).
        paths = []
        for set_name in self.sets_to_run:
            path = get_default_diags_path(set_name, run_type, False)
            paths.append(path)

        self.cfg_path = paths

        # Convert the .cfg file(s) to parameter objects.
        params = self.parser.get_cfg_parameters(
            files_to_open=paths, argparse_vals_only=False
        )

        # Update parameter objects using subclass with default values.
        params_final = self._convert_params_to_subclass(params)

        return params_final

    def _convert_params_to_subclass(
        self, parameters: List[CoreParameter]
    ) -> List[CoreParameter]:
        new_params: List[CoreParameter] = []

        # For each of the params, add in the default values using the parameter
        # classes in SET_TO_PARAMETERS.
        for param in parameters:
            set_key = param.sets[0]
            new_param = SET_TO_PARAMETERS[set_key]() + param

            new_params.append(new_param)

        return new_params

    def _get_parameters_with_api_only(
        self, parameters: List[CoreParameter]
    ) -> List[CoreParameter]:
        """Get the parameters defined through the Python API only.

        This method replaces CoreParameter objects with the related sub-class
        for the specified set.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of parameter objects.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects defined through the Python API only.
            Any non-CoreParameter objects will be replaced by a sub-class based
            on the set (``SETS_TO_PARAMETERS``).
        """
        params = []

        if len(self.sets_to_run) == 0:
            sets_to_run = [param.sets for param in parameters]
            self.sets_to_run = list(chain.from_iterable(sets_to_run))

        for set_name in self.sets_to_run:
            # For each of the set_names, get the corresponding parameter.
            param = self._get_instance_of_param_class(
                SET_TO_PARAMETERS[set_name], parameters
            )

            # Since each parameter will have lots of default values, we want to remove them.
            # Otherwise when calling get_parameters(), these default values
            # will take precedence over values defined in other_params.
            self._remove_attrs_with_default_values(param)
            param.sets = [set_name]

            # Makes sure that any parameters that are selectors will be in param.
            self._add_attrs_with_default_values(param)

            params.append(param)

        self.parser.check_values_of_params(params)

        return params

    def _add_parent_attrs_to_children(self, parameters):
        """
        For any parameter class that's inherited from another, copy
        the attributes of the parent to the child.

        Ex: If the user wants to run set-specific parameters for
        'zonal_mean_2d', they'd pass in a ZonalMean2dParameter
        and a CoreParameter.
        But most of the important parameters are in CoreParameter,
        so copy them over to the ZonalMean2dParameter.
        """

        def get_parent(param):
            """
            From parameters, get any object that's a parent
            type to param.

            Ex: CoreParameter object is a parent of AreaMeanTimeSeriesParameter object
            """
            try:
                parent_class = param.__class__.__mro__[1]
                parent = self._get_instance_of_param_class(parent_class, parameters)
            except RuntimeError:
                parent = None

            return parent

        for i in range(len(parameters)):
            parent = get_parent(parameters[i])
            # Make sure that the new object is actually a parent.
            if not parent or type(parent) == type(parameters[i]):  # noqa:E721
                continue

            # Otherwise, add the the parent's attributes.
            # Since we're modifying this parent object (by
            # removing the default values before addition)
            # make a deepcopy first.
            parent = copy.deepcopy(parent)
            # find attributes that are not defaults

            nondefault_param_parent = self._find_attrs_with_nondefault_values(parent)
            nondefault_param_child = self._find_attrs_with_nondefault_values(
                parameters[i]
            )

            self._remove_attrs_with_default_values(parent)

            # Simply copy over all attribute from parent to children
            # parameters[i] += parent

            for attr in dir(parent):
                if not attr.startswith("_") and not hasattr(parameters[i], attr):
                    # This attr of parent is a user-defined one and does not
                    # already exist in the parameters[i] parameter object.
                    attr_value = getattr(parent, attr)
                    setattr(parameters[i], attr, attr_value)

            logger.info(
                list(set(nondefault_param_parent) - set(nondefault_param_child))
            )
            for attr in list(
                set(nondefault_param_parent) - set(nondefault_param_child)
            ):
                # 'seasons' is a corner case that don't need to get in to none-core sets, Ex. area mean time series
                if attr != "seasons":
                    attr_value = getattr(parent, attr)
                    setattr(parameters[i], attr, attr_value)

    def _add_attrs_with_default_values(self, param):
        """
        In the param, add any missing parameters
        with their default value.
        """
        new_instance = param.__class__()
        for attr in dir(new_instance):
            # Ignore any of the hidden attributes.
            if attr.startswith("_"):
                continue

            if not hasattr(param, attr):
                val = getattr(new_instance, attr)
                setattr(param, attr, val)

    def _remove_attrs_with_default_values(self, param):
        """
        In the param, remove any parameters that
        have their default value.
        """
        new_instance = param.__class__()
        for attr in dir(param):
            # Ignore any of the hidden attributes.
            if attr.startswith("_"):
                continue

            if hasattr(new_instance, attr) and getattr(new_instance, attr) == getattr(
                param, attr
            ):
                delattr(param, attr)

    def _find_attrs_with_nondefault_values(self, param):
        """
        In the param, find any parameters that
        have nondefault value.
        """
        nondefault_attr = []
        new_instance = param.__class__()
        for attr in dir(param):
            # Ignore any of the hidden attributes.
            if attr.startswith("_"):
                continue

            if hasattr(new_instance, attr) and getattr(new_instance, attr) != getattr(
                param, attr
            ):  # This is only valid when the attr values are lists not numpy array
                nondefault_attr.append(attr)
        return nondefault_attr

    def _get_instance_of_param_class(self, cls, parameters):
        """
        In the list of parameters, get the class for
        the parameter object corresponding to cls.
        """
        # Get the Method Resolution Order (MRO) for this class.
        # So get the list of the classes in the inheritance ordering.

        # Ex: For the 'zonal_mean_2d' set, ZonalMean2dParameter is
        # the parameter for it.
        # But if a user doesn't want to modify the set-specific
        # parameters for 'zonal_mean_2d', they can just pass in
        # a single CoreParameter object to run_diags().
        # Using this, we handle this use-case.
        class_types = cls.__mro__

        for cls_type in class_types:
            for p in parameters:
                # NOTE: This conditional is used instead of
                # `isinstance(p, cls_type)` because we want to check for exact
                # type matching and exclude sub-class matching.
                if type(p) is cls_type:
                    return p

        msg = "There's weren't any class of types {} in your parameters."
        raise RuntimeError(msg.format(class_types))


runner = Run()
