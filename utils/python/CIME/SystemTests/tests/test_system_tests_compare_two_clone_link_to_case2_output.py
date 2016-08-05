#!/usr/bin/env python

"""
This module contains unit tests of the method
SystemTestsCompareTwoClone._link_to_case2_output
"""

import unittest
import os
import shutil
import tempfile
from CIME.SystemTests.system_tests_compare_two_clone import SystemTestsCompareTwoClone

class TestLinkToCase2Output(unittest.TestCase):

    # ========================================================================
    # Test helper functions
    # ========================================================================

    def setUp(self):
        self._run2_suffix = 'run2'
        self._casename1 = 'mytest'
        self._casename2 = 'mytest.run2'
        self._rundir1 = tempfile.mkdtemp()
        self._rundir2 = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._rundir1, ignore_errors=True)
        shutil.rmtree(self._rundir2, ignore_errors=True)


    # ========================================================================
    # Begin actual tests
    # ========================================================================

    def test_basic(self):
        # Setup
        filename1 = '%s.clm2.h0.nc.%s'%(self._casename2, self._run2_suffix)
        filepath1 = os.path.join(self._rundir2, filename1)
        open(filepath1, 'w').close()

        filename2 = '%s.clm2.h1.nc.%s'%(self._casename2, self._run2_suffix)
        filepath2 = os.path.join(self._rundir2, filename2)
        open(filepath2, 'w').close()

        # Exercise
        SystemTestsCompareTwoClone._link_to_case2_output(
            casename1 = self._casename1,
            casename2 = self._casename2,
            rundir1 = self._rundir1,
            rundir2 = self._rundir2,
            run2suffix = self._run2_suffix)

        # Verify
        expected_link_filename1 = '%s.clm2.h0.nc.%s'%(self._casename1, self._run2_suffix)
        expected_link_filepath1 = os.path.join(self._rundir1, expected_link_filename1)
        self.assertTrue(os.path.islink(expected_link_filepath1))
        self.assertEqual(filepath1, os.readlink(expected_link_filepath1))

        expected_link_filename2 = '%s.clm2.h1.nc.%s'%(self._casename1, self._run2_suffix)
        expected_link_filepath2 = os.path.join(self._rundir1, expected_link_filename2)
        self.assertTrue(os.path.islink(expected_link_filepath2))
        self.assertEqual(filepath2, os.readlink(expected_link_filepath2))

    def test_existing_link(self):
        # Setup
        filename1 = '%s.clm2.h0.nc.%s'%(self._casename2, self._run2_suffix)
        filepath1 = os.path.join(self._rundir2, filename1)
        open(filepath1, 'w').close()

        # Create initial link via a call to _link_to_case2_output
        SystemTestsCompareTwoClone._link_to_case2_output(
            casename1 = self._casename1,
            casename2 = self._casename2,
            rundir1 = self._rundir1,
            rundir2 = self._rundir2,
            run2suffix = self._run2_suffix)

        # Exercise
        # See what happens when we try to recreate that link
        SystemTestsCompareTwoClone._link_to_case2_output(
            casename1 = self._casename1,
            casename2 = self._casename2,
            rundir1 = self._rundir1,
            rundir2 = self._rundir2,
            run2suffix = self._run2_suffix)

        # (No verification: Test passes if no exception was raised)
