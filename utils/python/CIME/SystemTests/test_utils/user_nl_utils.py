"""
This module contains functions for working with user_nl files in system tests.

Typical usage for a test that needs to tweak a component's user_nl files for
each of multiple runs is:

In the pre-build phase:

user_nl_utils.save_user_nl_files(caseroot, component)

Then, in the initial setup for each run:

user_nl_utils.append_to_saved_files(caseroot, component, contents)
"""

import shutil
import os
import glob

_DEFAULT_SAVE_DIRNAME = "saved_user_nl_files"

def save_user_nl_files(caseroot, component,
                       save_dirname = _DEFAULT_SAVE_DIRNAME):
    """
    Save original user_nl files so that these originals can be restored later.

    If copies already exist in the given save directory, they will NOT be
    overwritten. Thus, it is safe to call this multiple times with the same
    arguments: only the first call will have an effect. (However, this behavior
    means that, in the unlikely event that user_nl files were present in the
    given save directory before the first call to save_user_nl_files - such as
    due to a poorly-chosen save_dirname - then these original files will be kept
    in place, which is likely not what you want.)

    Arguments:

    casedir: (string) Full path to the case directory

    component: (string) name of component (e.g., 'clm'). This is used to determine
    which user_nl files are copied and later modified. For example, for
    component='clm', this object will operate on all user_nl files matching
    the pattern 'user_nl_clm*'. (We do a wildcard match to handle
    multi-instance tests.)

    save_dirname: (string) name of directory to be created within the case directory,
    containing saved copies of the relevant user_nl file(s)
    """

    save_fullpath = os.path.join(caseroot, save_dirname)

    if not os.path.exists(save_fullpath):
        os.makedirs(save_fullpath)

    files = _get_list_of_user_nl_files(caseroot, component)

    if len(files) == 0:
        raise RuntimeError('No user_nl files found for component ' + component)

    for one_file in files:
        orig_file = os.path.join(caseroot, one_file)
        saved_file = os.path.join(save_fullpath, one_file)
        if not os.path.exists(saved_file):
            shutil.copy(orig_file, saved_file)

def append_to_saved_files(caseroot, component, contents,
                          save_dirname = _DEFAULT_SAVE_DIRNAME):
    """
    Copy the saved files back to the case directory, then append the
    string given by 'contents' to the end of each of the saved
    user_nl files for the given component (there may be multiple such user_nl files in
    the case of a multi-instance test).

    Also puts new lines before and after the appended text - so 'contents'
    does not need to contain a trailing new line (but it's also okay if it
    does).

    Because this method starts with the saved version of the files (from the
    call to save_user_nl_files), you can NOT use this twice to append
    additional text: the second call will overwrite the contents added in
    the first call.

    This should be called after save_user_nl_files has been called with the same
    caseroot, component and save_dirname arguments.

    Arguments:

    casedir: (string) Full path to the case directory

    component: (string) name of component (e.g., 'clm'). This is used to determine
    which user_nl files are copied and later modified. For example, for
    component='clm', this object will operate on all user_nl files matching
    the pattern 'user_nl_clm*'. (We do a wildcard match to handle
    multi-instance tests.)

    contents: (string) contents to append to the end of each user_nl file

    save_dirname: (string) name of sub-directory within the case directory that
    contains saved copies of the relevant user_nl file(s)
    """

    save_fullpath = os.path.join(caseroot, save_dirname)

    files = _get_list_of_user_nl_files(save_fullpath, component)

    if len(files) == 0:
        raise RuntimeError('No user_nl files found for component ' + component
                           + ' in ' + save_fullpath)

    for one_file in files:
        saved_file = os.path.join(save_fullpath, one_file)
        new_file = os.path.join(caseroot, one_file)
        shutil.copy(saved_file, new_file)
        with open(new_file, 'a') as user_nl_file:
            user_nl_file.write('\n' + contents + '\n')


def _get_list_of_user_nl_files(path, component):
    """Get a list of all user_nl files in the current path for the component
    of interest. For a component 'foo', we match all files of the form
    user_nl_foo* - with a wildcard match at the end in order to match files
    in a multi-instance case.

    The list of returned files gives just the file names themselves (i.e.,
    the basenames).
    """

    file_pattern = 'user_nl_' + component + '*'
    file_list = glob.glob(os.path.join(path, file_pattern))
    file_basename_list = [os.path.basename(one_file) for one_file in file_list]

    return file_basename_list
