#!/usr/bin/env python

import unittest
import os
import shutil
import tempfile
from CIME.SystemTests.test_utils.user_nl_copier import UserNLCopier

class TestUserNLCopier(unittest.TestCase):

    # ========================================================================
    # Test helper functions
    # ========================================================================

    def setUp(self):
        self._casedir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._casedir, ignore_errors=True)

    def write_user_nl_file(self, component, contents, suffix=''):
        """Write contents to a user_nl file in the case directory. Returns the
        basename (i.e., not the full path) of the file that is created.

        For a component foo, with the default suffix of '', the file name will
        be user_nl_foo

        If the suffix is '_0001', the file name will be user_nl_foo_0001
        """

        filename = 'user_nl_' + component + suffix

        with open(os.path.join(self._casedir, filename), 'w') as user_nl_file:
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
        nl_copier = UserNLCopier(casedir=self._casedir, component=component,
                                 save_dirname=save_dirname)

        # Exercise
        nl_copier.save_user_nl_files()

        # Verify
        self.assertTrue(os.path.isfile(
            os.path.join(self._casedir, save_dirname, filename)),
            msg = 'copied file should exist in save directory')

    def test_save_with_no_file_raises_exception(self):
        # Setup
        nl_copier = UserNLCopier(casedir=self._casedir, component='foo')

        # Exercise & verify
        self.assertRaisesRegexp(RuntimeError, "No user_nl files found",
                                nl_copier.save_user_nl_files)

    def test_save_called_twice_raises_exception(self):
        # Define some variables
        component = 'foo'
        orig_contents = 'bar = 42\n'

        # Setup
        with open(os.path.join(self._casedir, 'user_nl_' + component), 'w') as user_nl_file:
            user_nl_file.write(orig_contents)

        nl_copier = UserNLCopier(casedir=self._casedir, component=component)

        # Exercise & Verify
        nl_copier.save_user_nl_files()
        self.assertRaises(RuntimeError, 
                          nl_copier.save_user_nl_files)

    def test_save_with_existing_directory(self):
        # It should be okay to call save with an existing directory.

        # Define some variables
        component = 'foo'
        save_dirname = 'saved_files'
        orig_contents = 'bar = 42\n'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)
        os.makedirs(os.path.join(self._casedir, save_dirname))
        nl_copier = UserNLCopier(casedir=self._casedir, component=component,
                                 save_dirname=save_dirname)

        # Exercise
        nl_copier.save_user_nl_files()

        # Verify
        self.assertTrue(os.path.isfile(
            os.path.join(self._casedir, save_dirname, filename)),
            msg = 'copied file should exist in save directory')
        
    def test_save_overwriting_existing_files_raises_exception(self):
        # Trying to call save in a way that would overwrite existing files
        # should raise an exception

        # Define some variables
        component = 'foo'
        save_dirname = 'saved_files'
        orig_contents = 'bar = 42\n'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)
        os.makedirs(os.path.join(self._casedir, save_dirname))
        shutil.copy(os.path.join(self._casedir, filename),
                    os.path.join(self._casedir, save_dirname, filename))
        nl_copier = UserNLCopier(casedir=self._casedir, component=component,
                                 save_dirname=save_dirname)

        # Exercise & verify
        self.assertRaisesRegexp(RuntimeError,
                                "Attempt to overwrite existing saved files",
                                nl_copier.save_user_nl_files)

    def test_save_multiple_files(self):
        # Define some variables
        component = 'foo'
        save_dirname = 'saved_files'
        orig_contents1 = 'bar = 42\n'
        orig_contents2 = 'bar = 17\n'

        # Setup
        filename1 = self.write_user_nl_file(component, orig_contents1, suffix='_0001')
        filename2 = self.write_user_nl_file(component, orig_contents2, suffix='_0002')

        nl_copier = UserNLCopier(casedir=self._casedir, component=component,
                                 save_dirname=save_dirname)

        # Exercise
        nl_copier.save_user_nl_files()

        # Verify
        self.assertTrue(os.path.isfile(
            os.path.join(self._casedir, save_dirname, filename1)),
            msg = 'copied file 1 should exist in save directory')
        self.assertTrue(os.path.isfile(
            os.path.join(self._casedir, save_dirname, filename2)),
            msg = 'copied file 2 should exist in save directory')


    def test_append(self):
        # Define some variables
        component = 'foo'
        # deliberately exclude new line from file contents, to make sure that's
        # handled correctly
        orig_contents = 'bar = 42'
        contents_to_append = 'baz = 101'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)
        nl_copier = UserNLCopier(casedir=self._casedir, component=component)
        nl_copier.save_user_nl_files()

        # Exercise
        nl_copier.append_to_saved_files(contents_to_append)

        # Verify
        expected_contents = orig_contents + '\n' + contents_to_append + '\n'
        self.assertFileContentsEqual(expected_contents,
                                     os.path.join(self._casedir, filename))


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
        nl_copier = UserNLCopier(casedir=self._casedir, component=component)
        nl_copier.save_user_nl_files()

        # First append should not affect final result
        nl_copier.append_to_saved_files(contents_to_append_first)

        # Exercise
        nl_copier.append_to_saved_files(contents_to_append_second)

        # Verify
        expected_contents1 = orig_contents1 + '\n' + contents_to_append_second + '\n'
        expected_contents2 = orig_contents2 + '\n' + contents_to_append_second + '\n'
        self.assertFileContentsEqual(expected_contents1,
                                     os.path.join(self._casedir, filename1))
        self.assertFileContentsEqual(expected_contents2,
                                     os.path.join(self._casedir, filename2))

    def test_append_without_save_raises_exception(self):
        # Setup
        nl_copier = UserNLCopier(casedir=self._casedir, component='foo')

        # Exercise & verify
        self.assertRaisesRegexp(RuntimeError,
                                "save_user_nl_files must be called before append_to_saved_files",
                                nl_copier.append_to_saved_files, 'bar')
