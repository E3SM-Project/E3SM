"""
This module contains functions for working with user_nl files in system tests.
"""

import os
import glob

def append_to_user_nl_files(caseroot, component, contents):
    """
    Append the string given by 'contents' to the end of each user_nl file for
    the given component (there may be multiple such user_nl files in the case of
    a multi-instance test).

    Also puts new lines before and after the appended text - so 'contents'
    does not need to contain a trailing new line (but it's also okay if it
    does).

    Args:
        caseroot (str): Full path to the case directory

        component (str): Name of component (e.g., 'clm'). This is used to
            determine which user_nl files are appended to. For example, for
            component='clm', this function will operate on all user_nl files
            matching the pattern 'user_nl_clm*'. (We do a wildcard match to
            handle multi-instance tests.)

        contents (str): Contents to append to the end of each user_nl file
    """

    files = _get_list_of_user_nl_files(caseroot, component)

    if len(files) == 0:
        raise RuntimeError('No user_nl files found for component ' + component)

    for one_file in files:
        with open(one_file, 'a') as user_nl_file:
            user_nl_file.write('\n' + contents + '\n')

def _get_list_of_user_nl_files(path, component):
    """Get a list of all user_nl files in the current path for the component
    of interest. For a component 'foo', we match all files of the form
    user_nl_foo* - with a wildcard match at the end in order to match files
    in a multi-instance case.

    The list of returned files gives their full path.
    """

    file_pattern = 'user_nl_' + component + '*'
    file_list = glob.glob(os.path.join(path, file_pattern))

    return file_list
