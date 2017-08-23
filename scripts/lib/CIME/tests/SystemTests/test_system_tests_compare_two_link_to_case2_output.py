#!/usr/bin/env python

"""
This module contains unit tests of the method
SystemTestsCompareTwo._link_to_case2_output
"""

# Ignore privacy concerns for unit tests, so that unit tests can access
# protected members of the system under test
#
# pylint:disable=protected-access

import unittest
import os
import shutil
import tempfile
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.tests.case_fake import CaseFake

# ========================================================================
# Fake version of SystemTestsCompareTwo that overrides some functionality for
# the sake of unit testing
# ========================================================================

class SystemTestsCompareTwoFake(SystemTestsCompareTwo):
    def __init__(self,
                 case1,
                 run_two_suffix = 'test'):

        SystemTestsCompareTwo.__init__(
            self,
            case1,
            separate_builds = False,
            run_two_suffix = run_two_suffix)

    # ------------------------------------------------------------------------
    # Stubs of methods called by SystemTestsCommon.__init__ that interact with
    # the system or case object in ways we want to avoid here
    # ------------------------------------------------------------------------

    def _init_environment(self, caseroot):
        pass

    def _init_locked_files(self, caseroot, expected):
        pass

    def _init_case_setup(self):
        pass

    # ------------------------------------------------------------------------
    # Stubs of methods that are typically provided by the individual test
    # ------------------------------------------------------------------------

    def _case_one_setup(self):
        pass

    def _case_two_setup(self):
        pass

# ========================================================================
# Test class itself
# ========================================================================

class TestLinkToCase2Output(unittest.TestCase):

    # ========================================================================
    # Test helper functions
    # ========================================================================

    def setUp(self):
        self.original_wd = os.getcwd()
        # Create a sandbox in which case directories can be created
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        # Some tests trigger a chdir call in the SUT; make sure we return to the
        # original directory at the end of the test
        os.chdir(self.original_wd)

        shutil.rmtree(self.tempdir, ignore_errors=True)

    def setup_test_and_directories(self, casename1, run2_suffix):
        """
        Returns test object
        """

        case1root = os.path.join(self.tempdir, casename1)
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1, run_two_suffix = run2_suffix)
        mytest._case1.make_rundir()  #pylint: disable=maybe-no-member
        mytest._case2.make_rundir()  #pylint: disable=maybe-no-member

        return mytest

    def create_file_in_rundir2(self, mytest, core_filename, run2_suffix):
        """
        Creates a file in rundir2 named CASE2.CORE_FILENAME.nc.RUN2_SUFFIX
        (where CASE2 is the casename of case2)

        Returns full path to the file created
        """
        filename = '%s.%s.nc.%s'%(
            mytest._case2.get_value('CASE'),
            core_filename,
            run2_suffix)
        filepath = os.path.join(mytest._case2.get_value('RUNDIR'), filename)
        open(filepath, 'w').close()
        return filepath

    # ========================================================================
    # Begin actual tests
    # ========================================================================

    def test_basic(self):
        # Setup
        casename1 = 'mytest'
        run2_suffix = 'run2'

        mytest = self.setup_test_and_directories(casename1, run2_suffix)
        filepath1 = self.create_file_in_rundir2(mytest, 'clm2.h0', run2_suffix)
        filepath2 = self.create_file_in_rundir2(mytest, 'clm2.h1', run2_suffix)

        # Exercise
        mytest._link_to_case2_output()

        # Verify
        expected_link_filename1 = '%s.clm2.h0.nc.%s'%(casename1, run2_suffix)
        expected_link_filepath1 = os.path.join(mytest._case1.get_value('RUNDIR'),
                                               expected_link_filename1)
        self.assertTrue(os.path.islink(expected_link_filepath1))
        self.assertEqual(filepath1, os.readlink(expected_link_filepath1))

        expected_link_filename2 = '%s.clm2.h1.nc.%s'%(casename1, run2_suffix)
        expected_link_filepath2 = os.path.join(mytest._case1.get_value('RUNDIR'),
                                               expected_link_filename2)
        self.assertTrue(os.path.islink(expected_link_filepath2))
        self.assertEqual(filepath2, os.readlink(expected_link_filepath2))

    def test_existing_link(self):
        # Setup
        casename1 = 'mytest'
        run2_suffix = 'run2'

        mytest = self.setup_test_and_directories(casename1, run2_suffix)
        self.create_file_in_rundir2(mytest, 'clm2.h0', run2_suffix)

        # Create initial link via a call to _link_to_case2_output
        mytest._link_to_case2_output()

        # Exercise
        # See what happens when we try to recreate that link
        mytest._link_to_case2_output()

        # (No verification: Test passes if no exception was raised)
