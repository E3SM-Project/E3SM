import os
import cdp.cdp_run
import acme_diags
from acme_diags.acme_diags_driver import main
from acme_diags.parser.core_parser import CoreParser
from acme_diags.parameter import SET_TO_PARAMETERS
from acme_diags.parameter.core_parameter import CoreParameter


class Run():
    """
    Used to run diagnostics.
    A class is needed because we often need to store some
    state regarding what we need to run, like the sets selected.
    """
    def __init__(self):
        self.sets_to_run = CoreParameter().sets
        self.parser = CoreParser()

    def run_diags(self, parameters):
        """
        Based on sets_to_run, run the diags with the list of parameters.
        """
        if not parameters or not isinstance(parameters, list):
            msg = 'You must pass in a list of parameter objects.'
            raise RuntimeError(msg)
        
        # For each of the passed in parameters, we can only have one of
        # each type.
        types = set([p.__class__ for p in parameters])
        if len(types) != len(parameters):
            msg = 'You passed in two or more parameters of the same type.'
            raise RuntimeError(msg)

        final_params = []
        other_params = self._get_other_diags(parameters[0].run_type)

        for set_name in self.sets_to_run:
            # For each of the set_names, corresponding parameter.
            param = self._get_instance_of_param_class(set_name, parameters)
            self._remove_attrs_with_default_values(param)
            param.sets = [set_name]

            params = self.parser.get_parameters(orig_parameters=param, other_parameters=other_params,
                cmd_default_vars=False, argparse_vals_only=False)
            # The select() call in get_parameters() was made for the original
            # command-line way of using CDP.
            # We just call it manually with the parameter object param.
            params = self.parser.select(param, params)
            final_params.extend(params)

        self.parser.check_values_of_params(final_params)

        main(final_params)


    def _remove_attrs_with_default_values(self, param):
        """
        In the param, remove any parameters that have their default value.
        However, since `selectors` and `granulate` are special parameters
        (from these, more parameter objects are created), don't remove them.
        """
        new_instance = param.__class__()
        for attr in dir(param):
            # Ignore any of the hidden attributes.
            if attr.startswith('_') or attr in ['selectors', 'granulate']:
                continue
            
            if hasattr(new_instance, attr) and \
                getattr(new_instance, attr) == getattr(param, attr):
                delattr(param, attr)


    def _get_instance_of_param_class(self, set_name, parameters):
        """
        In the list of parameters, get the class for
        the parameter object corresponding to set_name.
        """
        class_type = SET_TO_PARAMETERS[set_name]

        for p in parameters:
            if isinstance(p, class_type):
                return p
        
        msg = "There's no class of type {} in your parameters."
        raise RuntimeError(msg.format(class_type))


    def _get_default_diags_path(self, set_name, run_type):
        """
        Returns the path for the default diags for plotset set_name.
        These are different depending on the run_type.
        """
        folder = '{}'.format(set_name)
        fnm = '{}_{}.cfg'.format(set_name, run_type)
        pth = os.path.join(acme_diags.INSTALL_PATH, folder, fnm)

        print('Using {} for {}.'.format(pth, set_name))
        if not os.path.exists(pth):
            raise RuntimeError(
                "Plotting via set '{}' not supported, file {} not installed".format(set_name, fnm))
        return pth


    def _get_other_diags(self, run_type):
        """
        If the user has ran the script with a -d, get the diags for that.
        If not, load the default diags based on sets_to_run and run_type.
        """
        args = self.parser.view_args()

        # If the user has passed in args with -d.
        if args.other_parameters:
            return self.parser.get_other_parameters(argparse_vals_only=False)
        else:
            default_diags_paths = [self._get_default_diags_path(set_name, run_type) for set_name in self.sets_to_run]
            return self.parser.get_other_parameters(files_to_open=default_diags_paths, argparse_vals_only=False)
