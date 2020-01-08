import os
import copy
import cdp.cdp_run
import acme_diags
from acme_diags.acme_diags_driver import main, get_default_diags_path
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
        final_params = self.get_final_parameters(parameters)
        if not final_params:
            msg = 'No parameters we able to be extracted.'
            msg += ' Please check the parameters you defined.'
            raise RuntimeError(msg)

        main(final_params)


    def get_final_parameters(self, parameters):
        """
        Based on sets_to_run and the list of parameters,
        get the final list of paremeters to run the diags on.
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

        self._add_parent_attrs_to_children(parameters)

        final_params = []

        for set_name in self.sets_to_run:
            other_params = self._get_other_diags(parameters[0].run_type)

            # For each of the set_names, get the corresponding parameter.
            param = self._get_instance_of_param_class(SET_TO_PARAMETERS[set_name], parameters)
            

            # Since each parameter will have lots of default values, we want to remove them.
            # Otherwise when calling get_parameters(), these default values
            # will take precedence over values defined in other_params.
            self._remove_attrs_with_default_values(param)
            param.sets = [set_name]

            params = self.parser.get_parameters(orig_parameters=param, other_parameters=other_params,
                cmd_default_vars=False, argparse_vals_only=False)
           
            # Makes sure that any parameters that are selectors
            # will be in param.
            self._add_attrs_with_default_values(param)
            # The select() call in get_parameters() was made for the original
            # command-line way of using CDP.
            # We just call it manually with the parameter object param.
            params = self.parser.select(param, params)


            final_params.extend(params)

        self.parser.check_values_of_params(final_params)

        return final_params


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
            if not parent or type(parent) == type(parameters[i]):
                continue

            # Otherwise, add the the parent's attributes.
            # Since we're modifying this parent object (by
            # removing the default values before addition)
            # make a deepcopy first.
            parent = copy.deepcopy(parent)
            #find attributes that are not defaults

            nondefault_param_parent = self._find_attrs_with_nondefault_values(parent)
            nondefault_param_child = self._find_attrs_with_nondefault_values(parameters[i])


            self._remove_attrs_with_default_values(parent)
 

            #Simply copy over all attribute from parent to children
            #parameters[i] += parent
             
            for attr in dir(parent):
                if not attr.startswith('_') and not hasattr(parameters[i], attr):
                    # This attr of parent is a user-defined one and does not
                    # already exist in the parameters[i] parameter object.
                    attr_value = getattr(parent, attr)
                    setattr(parameters[i], attr, attr_value)

            print(list(set(nondefault_param_parent) - \
                set(nondefault_param_child)))
            for attr in list(set(nondefault_param_parent) - \
                set(nondefault_param_child)):
                #'seasons' is a corner case that don't need to get in to none-core sets, Ex. area mean time series
                if attr != 'seasons':
                    attr_value = getattr(parent, attr)
                    setattr(parameters[i], attr, attr_value)

        #for i in range(len(parameters)):
        #    attrs = vars(parameters[i])
        #    print('all parameters', ','.join("%s: %s" % item for item in attrs.items()))


    def _add_attrs_with_default_values(self, param):
        """
        In the param, add any missing parameters
        with their default value.
        """
        new_instance = param.__class__()
        for attr in dir(new_instance):
            # Ignore any of the hidden attributes.
            if attr.startswith('_'):
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
            if attr.startswith('_'):
                continue
            
            if hasattr(new_instance, attr) and \
                getattr(new_instance, attr) == getattr(param, attr):
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
            if attr.startswith('_'):
                continue
            
            if hasattr(new_instance, attr) and \
                getattr(new_instance, attr) != getattr(param, attr):   # This is only valid when the attr values are lists not numpy array
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
                if type(p) == cls_type:
                    return p
        
        msg = "There's weren\'t any class of types {} in your parameters."
        raise RuntimeError(msg.format(class_types))


    def _get_other_diags(self, run_type):
        """
        If the user has ran the script with a -d, get the diags for that.
        If not, load the default diags based on sets_to_run and run_type.
        """
        args = self.parser.view_args()

        # If the user has passed in args with -d.
        if args.other_parameters:
            params = self.parser.get_other_parameters(argparse_vals_only=False)
        else:
            default_diags_paths = [get_default_diags_path(set_name, run_type, False) for set_name in self.sets_to_run]
            params = self.parser.get_other_parameters(files_to_open=default_diags_paths, argparse_vals_only=False)

        # For each of the params, add in the default values
        # using the parameter classes in SET_TO_PARAMETERS.
        for i in range(len(params)):
            params[i] = SET_TO_PARAMETERS[params[i].sets[0]]() + params[i]
        
        return params


runner = Run()

