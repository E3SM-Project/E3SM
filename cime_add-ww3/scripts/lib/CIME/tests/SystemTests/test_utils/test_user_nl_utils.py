#!/usr/bin/env python

import unittest
import os
import shutil
import tempfile
from CIME.SystemTests.test_utils import user_nl_utils
import six

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

    def test_append(self):
        # Define some variables
        component = 'foo'
        # deliberately exclude new line from file contents, to make sure that's
        # handled correctly
        orig_contents = 'bar = 42'
        contents_to_append = 'baz = 101'

        # Setup
        filename = self.write_user_nl_file(component, orig_contents)

        # Exercise
        user_nl_utils.append_to_user_nl_files(caseroot = self._caseroot,
                                              component = component,
                                              contents = contents_to_append)

        # Verify
        expected_contents = orig_contents + '\n' + contents_to_append + '\n'
        self.assertFileContentsEqual(expected_contents,
                                     os.path.join(self._caseroot, filename))

    def test_append_multiple_files(self):
        # Simulates a multi-instance test
        component = 'foo'
        orig_contents1 = 'bar = 42'
        orig_contents2 = 'bar = 17'
        contents_to_append = 'baz = 101'

        # Setup
        filename1 = self.write_user_nl_file(component, orig_contents1, suffix='_0001')
        filename2 = self.write_user_nl_file(component, orig_contents2, suffix='_0002')

        # Exercise
        user_nl_utils.append_to_user_nl_files(caseroot = self._caseroot,
                                              component = component,
                                              contents = contents_to_append)

        # Verify
        expected_contents1 = orig_contents1 + '\n' + contents_to_append + '\n'
        expected_contents2 = orig_contents2 + '\n' + contents_to_append + '\n'
        self.assertFileContentsEqual(expected_contents1,
                                     os.path.join(self._caseroot, filename1))
        self.assertFileContentsEqual(expected_contents2,
                                     os.path.join(self._caseroot, filename2))


    def test_append_without_files_raises_exception(self):
        # This test verifies that you get an exception if you call
        # append_to_user_nl_files when there are no user_nl files of interest

        # Define some variables
        component_exists = 'foo'
        component_for_append = 'bar'

        # Setup
        # Create file in caseroot for component_exists, but not for component_for_append
        self.write_user_nl_file(component_exists, 'irrelevant contents')

        # Exercise & verify
        six.assertRaisesRegex(self, RuntimeError, "No user_nl files found",
                                user_nl_utils.append_to_user_nl_files,
                                caseroot = self._caseroot,
                                component = component_for_append,
                                contents = 'irrelevant contents to append')
