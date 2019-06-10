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
                
        other_parameters = self._get_other_diags(parameters[0].run_type)

        # print([p.sets for p in other_parameters])
        # print([p.variables for p in other_parameters])

        for set_name in self.sets_to_run:
            # For each of the set_names, corresponding parameter.
            param_class = SET_TO_PARAMETERS[set_name]
            param = self._get_instance_of_param_class(param_class, parameters)

            # Combine the key-value pairs from param into other_parameters.
            self._combine_params(param, other_parameters)

        # for _ in range(5):
        #     print('*'*40)
        # print([p.sets for p in other_parameters])
        # print([p.variables for p in other_parameters])

        main(other_parameters)


    def _remove_attrs_with_default_values(self, param):
        """
        In the param, remove any parameters that
        have their default value.
        """
        new_instance = param.__class__()
        for attr in dir(param):
            # Ignore any of the hidden attributes.
            if attr.startswith('_'):
                continue
            
            if hasattr(new_instance, attr) and \
                getattr(new_instance, attr) == getattr(param, attr):
                delattr(param, attr)


    def _add_attrs_with_default_values(self, param):
        """
        In the param, add in any missing attributes.
        When adding them, the original attribute is added.
        """
        new_instance = param.__class__()
        for attr in dir(new_instance):
            # Ignore any of the hidden attributes.
            if attr.startswith('_'):
                continue
            
            if not hasattr(param, attr):
                val = getattr(new_instance, attr)
                setattr(param, attr, val)


    def _combine_params(self, param, other_params):
        """
        Combine a single parameter object param with a list
        of parameter objects other_params.
        """
        cmdline_parameters = self.parser.get_cmdline_parameters(argparse_vals_only=False)
        
        # We don't want to add the selectors to each of the parameters.
        # Because if we do, it'll select all of the parameters at the end during the selection step.
        vars_to_ignore = self.parser._get_selectors(cmdline_parameters, param, other_params)
        # We need to remove any parameters that still have their default value.
        # Ex: param.variables == [] is the default value for the variables parameter.
        #     We need to remove this, otherwise after combining the parameters,
        #     each will have a value of variables == [].
        #     This doesn't make any sense.
        self._remove_attrs_with_default_values(param)
        print(vars(other_params[0]))
        self.parser.combine_params(orig_parameters=param, other_parameters=other_params, vars_to_ignore=vars_to_ignore)
        for param in other_params:
            # After combining, these default parameters need to be added back in.
            self._add_attrs_with_default_values(param)

        print(vars(other_params[0]))


    def _get_instance_of_param_class(self, class_type, parameters):
        """
        In the list of parameters, get the first class of instance class_type.
        """
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
