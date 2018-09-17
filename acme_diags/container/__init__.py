"""
Functionality used for when running e3sm_diags as a container.
"""
import os


def is_container():
    """
    Returns True if e3sm_diags is running as a container.
    """
    # If e3sm_diags has a current wording directory as the container's,
    # then it's probably running as a container.
    # Why would someone otherwise actually run the diags in a folder 'e3sm_diags_container_cwd'?
    return os.getcwd() == '/e3sm_diags_container_cwd'

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
