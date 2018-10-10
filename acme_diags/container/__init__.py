"""
Functionality used for when running e3sm_diags as a container.
"""
import os


def is_container():
    """
    Returns True if e3sm_diags is running as a container.
    """
    # Check if the environmental variable set in the Dockerfile
    # is present and is True.
    var = 'E3SM_DIAGS_CONTAINER'
    return var in os.environ and os.environ[var] == 'true'

def containerize_parameter(params_obj):
    """
    Given a set of parameters, modify them to be prepared
    to be ran in e3sm_diags as a container.
    """
    params = ['reference_data_path', 'test_data_path', 'results_dir']

    # When running in a container, set each parameter to '/{parameter}'.
    # This makes it easy for e3sm_diags as a container to get the input data
    # and output the results.
    # But make a backup of the values, so when saving the command in the viewer,
    # the results can be recreated with or without running a container.
    for p in params:
        orig_attr_value = getattr(params_obj, p, '')
        setattr(params_obj, 'orig_{}'.format(p), orig_attr_value)
        setattr(params_obj, p, '/{}'.format(p))

def decontainerize_parameter(params_obj):
    """
    Given a set of parameters already ran through containerize_parameter(),
    reverse what was done.
    """
    params = ['reference_data_path', 'test_data_path', 'results_dir']

    # Set each of the params back to their original value.
    for p in params:
        orig_attr_value = getattr(params_obj, 'orig_{}'.format(p))
        setattr(params_obj, p, orig_attr_value)
        delattr(params_obj, 'orig_{}'.format(p))
