import shutil
import os
import glob

class UserNLCopier(object):
    """This class can be used to save the original user_nl files for a given
    component, then restore them with some appended modifications for each run
    of a test.

    Usage is:

    In the constructor for a test:
    self._nl_copier = UserNLCopier(casedir, component)

    Then, after the user_nl files have been copied into the case directory:
    self._nl_copier.save_user_nl_files()

    Then, in the initial setup for each run:
    self._nl_copier.append_to_saved_files(string_to_append)


    If you need to handle user_nl files for multiple components in your test, you should
    create a separate UserNLCopier object for each component (e.g., one for clm,
    one for cam, etc.)
    """

    def __init__(self, casedir, component, save_dirname = 'user_nl_orig'):
        """Creates a UserNLCopier object

        casedir: full path of the case directory

        component: name of component (e.g., 'clm'). This is used to determine
        which user_nl files are copied and later modified. For example, for
        component='clm', this object will operate on all user_nl files matching
        the pattern 'user_nl_clm*'. (We do a wildcard match to handle
        multi-instance tests.)

        save_dirname: name of directory to be created within the case directory,
        containing saved copies of the relevant user_nl file(s)
        """

        self._casedir = casedir
        self._component = component
        self._save_dirname = save_dirname
        self._save_fullpath = os.path.join(self._casedir, self._save_dirname)
        self._files_saved = False

    def save_user_nl_files(self):
        """Save original user_nl files so that these originals can be restored
        later. This should be called exactly once per test, after the user_nl
        files have been copied into the case directory."""

        if self._files_saved:
            raise RuntimeError('Attempt to call save_user_nl_files twice')

        if not os.path.exists(self._save_fullpath):
            os.makedirs(self._save_fullpath)

        files = self._get_list_of_user_nl_files(self._casedir)

        if len(files) == 0:
            raise RuntimeError('No user_nl files found for component ' + self._component)

        for one_file in files:
            orig_file = os.path.join(self._casedir, one_file)
            saved_file = os.path.join(self._save_fullpath, one_file)
            if (os.path.exists(saved_file)):
                raise RuntimeError('Attempt to overwrite existing saved files in ' +
                                   self._save_fullpath)
            shutil.copy(orig_file, saved_file)

        self._files_saved = True

    def append_to_saved_files(self, contents):
        """Copy the saved files back to the case directory, then append the
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

        This must be called after save_user_nl_files.

        contents: string giving the contents to append to the end of each user_nl file
        """

        if not self._files_saved:
            raise RuntimeError('save_user_nl_files must be called before append_to_saved_files')

        files = self._get_list_of_user_nl_files(self._save_fullpath)

        if len(files) == 0:
            raise RuntimeError('No user_nl files found for component ' + self._component)

        for one_file in files:
            saved_file = os.path.join(self._save_fullpath, one_file)
            new_file = os.path.join(self._casedir, one_file)
            shutil.copy(saved_file, new_file)
            with open(new_file, 'a') as user_nl_file:
                user_nl_file.write('\n' + contents + '\n')


    def _get_list_of_user_nl_files(self, path):
        """Get a list of all user_nl files in the current path for the component
        of interest. For a component 'foo', we match all files of the form
        user_nl_foo* - with a wildcard match at the end in order to match files
        in a multi-instance case.

        The list of returned files gives just the file names themselves (i.e.,
        the basenames).
        """

        file_pattern = 'user_nl_' + self._component + '*'
        file_list = glob.glob(os.path.join(self._casedir, file_pattern))
        file_basename_list = [os.path.basename(one_file) for one_file in file_list]

        return file_basename_list
