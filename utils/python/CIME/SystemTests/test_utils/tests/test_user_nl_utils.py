#!/usr/bin/env python

import unittest
import os
import shutil
import tempfile
from CIME.SystemTests.test_utils import user_nl_utils

class TestUserNLCopier(unittest.TestCase):

    # ========================================================================
    # Test helper functions
    # ========================================================================

    def setUp(self):
        self._caseroot = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._caseroot, ignore_errors=True)

    def write_user_nl_file(self, component, contents, suffix=''):
        """Write contents to a user_nl file in the case directory. Returns the
        basename (i.e., not the full path) of the file that is created.

        For a component foo, with the default suffix of '', the file name will
        be user_nl_foo

        If the suffix is '_0001', the file name will be user_nl_foo_0001
        """

        filename = 'user_nl_' + component + suffix

        with open(os.path.join(self._caseroot, filename), 'w') as user_nl_file:
            user_nl_file.write(contents)

        return filename

    def assertFileContentsEqual(self, expected, filepath, msg=None):
        """Asserts that the contents of the file given by 'filepath' are equal to
        the string given by 'expected'. 'msg' gives an optional message to be
        printed if the assertion fails."""

        with open(filepath, 'r') as myfile:
            contents = myfile.read()

        self.assertEqual(expected, contents, msg=msg)

    # ========================================================================
    # Begin actual tests
    # ========================================================================

    def test_save(self):
        # Define some variables
        component = 'foo'
        save_dirname = 'saved_files'
        orig_contents = 'bar = 42\n'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)

        # Exercise
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component,
                                         save_dirname = save_dirname)

        # Verify
        self.assertTrue(os.path.isfile(
            os.path.join(self._caseroot, save_dirname, filename)),
            msg = 'copied file should exist in save directory')

    def test_save_multiple_files(self):
        # Define some variables
        component = 'foo'
        save_dirname = 'saved_files'
        orig_contents1 = 'bar = 42\n'
        orig_contents2 = 'bar = 17\n'

        # Setup
        filename1 = self.write_user_nl_file(component, orig_contents1, suffix='_0001')
        filename2 = self.write_user_nl_file(component, orig_contents2, suffix='_0002')

        # Exercise
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component,
                                         save_dirname = save_dirname)

        # Verify
        self.assertTrue(os.path.isfile(
            os.path.join(self._caseroot, save_dirname, filename1)),
            msg = 'copied file 1 should exist in save directory')
        self.assertTrue(os.path.isfile(
            os.path.join(self._caseroot, save_dirname, filename2)),
            msg = 'copied file 2 should exist in save directory')


    def test_save_with_existing_directory(self):
        # It should be okay to call save with an existing directory.

        # Define some variables
        component = 'foo'
        save_dirname = 'saved_files'
        orig_contents = 'bar = 42\n'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)
        os.makedirs(os.path.join(self._caseroot, save_dirname))

        # Exercise
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component,
                                         save_dirname = save_dirname)

        # Verify
        self.assertTrue(os.path.isfile(
            os.path.join(self._caseroot, save_dirname, filename)),
            msg = 'copied file should exist in save directory')

    def test_save_twice_should_not_overwrite_existing_files(self):
        # If you call save_user_nl_files when there are already saved user_nl
        # files in the save directory, then they should not be overwritten

        # Define some variables
        component = 'foo'
        save_dirname = 'saved_files'
        orig_contents = 'bar = 42\n'
        new_contents = 'bar = 17\n'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)
        # do an initial copy to the save directory to set things up for the real
        # test:
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component,
                                         save_dirname = save_dirname)
        # now overwrite the file in the case directory:
        filename = self.write_user_nl_file(component, new_contents)

        # Exercise
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component,
                                         save_dirname = save_dirname)

        # Verify
        self.assertFileContentsEqual(orig_contents,
                                     os.path.join(self._caseroot, save_dirname, filename))

    def test_save_with_no_file_raises_exception(self):
        # Exercise & verify
        self.assertRaisesRegexp(RuntimeError, "No user_nl files found",
                                user_nl_utils.save_user_nl_files,
                                caseroot = self._caseroot,
                                component = 'foo')


    def test_append(self):
        # Define some variables
        component = 'foo'
        # deliberately exclude new line from file contents, to make sure that's
        # handled correctly
        orig_contents = 'bar = 42'
        contents_to_append = 'baz = 101'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component)

        # Exercise
        user_nl_utils.append_to_saved_files(caseroot = self._caseroot,
                                            component = component,
                                            contents = contents_to_append)

        # Verify
        expected_contents = orig_contents + '\n' + contents_to_append + '\n'
        self.assertFileContentsEqual(expected_contents,
                                     os.path.join(self._caseroot, filename))


    def test_two_appends(self):
        # If you call append twice, only the second append should be present in
        # the final file.
        #
        # This test also tests appending with multi-instance

        # Define some variables
        component = 'foo'
        orig_contents1 = 'bar = 42'
        orig_contents2 = 'bar = 17'
        contents_to_append_first = 'baz = 101'
        contents_to_append_second = 'baz = 201'

        # Setup
        filename1 = self.write_user_nl_file(component, orig_contents1, suffix='_0001')
        filename2 = self.write_user_nl_file(component, orig_contents2, suffix='_0002')
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component)

        # First append should not affect final result
        user_nl_utils.append_to_saved_files(caseroot = self._caseroot,
                                            component = component,
                                            contents = contents_to_append_first)

        # Exercise
        user_nl_utils.append_to_saved_files(caseroot = self._caseroot,
                                            component = component,
                                            contents = contents_to_append_second)

        # Verify
        expected_contents1 = orig_contents1 + '\n' + contents_to_append_second + '\n'
        expected_contents2 = orig_contents2 + '\n' + contents_to_append_second + '\n'
        self.assertFileContentsEqual(expected_contents1,
                                     os.path.join(self._caseroot, filename1))
        self.assertFileContentsEqual(expected_contents2,
                                     os.path.join(self._caseroot, filename2))

    def test_append_without_savedir_raises_exception(self):
        # This test verifies that you get an exception if you try to call
        # append_to_saved_files without first calling save_user_nl_files

        # Exercise and verify
        self.assertRaisesRegexp(RuntimeError, "No user_nl files found",
                                user_nl_utils.append_to_saved_files,
                                caseroot = self._caseroot,
                                component = 'foo',
                                contents = 'bar')

    def test_append_without_saved_files_raises_exception(self):
        # This test verifies that you get an exception if you call
        # append_to_saved_files using a save directory that does not contain the
        # user_nl files of interest

        # Define some variables
        component_saved = 'foo'
        component_for_append = 'bar'

        # Setup
        # Create files in caseroot for both component_saved and component_for_append...
        filename = self.write_user_nl_file(component_saved, 'irrelevant contents')
        filename = self.write_user_nl_file(component_for_append, 'other irrelevant contents')
        # ... but only call save on component_saved
        user_nl_utils.save_user_nl_files(caseroot = self._caseroot,
                                         component = component_saved)

        # Exercise & verify
        self.assertRaisesRegexp(RuntimeError, "No user_nl files found",
                                user_nl_utils.append_to_saved_files,
                                caseroot = self._caseroot,
                                component = component_for_append,
                                contents = 'irrelevant contents to append')

