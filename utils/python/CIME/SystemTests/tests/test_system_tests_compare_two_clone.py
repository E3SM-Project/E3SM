#!/usr/bin/env python

"""
This module contains unit tests of the core logic in SystemTestsCompareTwoClone.
"""

import unittest
from CIME.SystemTests.system_tests_compare_two_clone import SystemTestsCompareTwoClone

import os
import shutil
import tempfile
from copy import deepcopy

# ========================================================================
# Fake version of Case object that provides the functionality needed for these
# tests
#
# TODO(wjs, 2016-08-10) Should this be moved into case.py? (Is it useful enough
# to warrant that?)
# ========================================================================

class CaseFake(object):
    def __init__(self, case_root, create_case_root=True):
        """
        Initialize a new case object for the given case_root directory.

        Args:
            case_root (str): path to CASEROOT
            create_case_root (bool): If True, creates the directory given by case_root
        """
        self.vars = dict()
        if create_case_root:
            os.makedirs(case_root)
        self.set_value('CASEROOT', case_root)
        casename = os.path.basename(case_root)
        self.set_value('CASE', casename)
        self.set_value('CASEBASEID', casename)

    def get_value(self, item):
        """
        Get the value of the given item

        Returns None if item isn't set for this case

        Args:
            item (str): variable of interest
        """
        return self.vars.get(item)

    def set_value(self, item, value):
        """
        Set the value of the given item to the given value

        Args:
            item (str): variable of interest
            value (any type): new value for item
        """
        self.vars[item] = value

    def copy(self, newcasename, newcaseroot):
        """
        Create and return a copy of self, but with CASE and CASEBASEID set to newcasename
        and CASEROOT set to newcaseroot

        Args:
            newcasename (str): new value for CASE
            newcaseroot (str): new value for CASEROOT
        """
        newcase = deepcopy(self)
        newcase.set_value('CASE', newcasename)
        newcase.set_value('CASEBASEID', newcasename)
        newcase.set_value('CASEROOT', newcaseroot)

        return newcase

    def create_clone(self, newcase, keepexe=False):
        """
        Create a clone of the current case. Also creates the CASEROOT directory
        for the clone case (given by newcase).

        Args:
            newcase (str): full path to the new case. This directory should not
                already exist; it will be created
            keepexe (bool, optional): Ignored

        Returns the clone case object
        """
        newcaseroot = os.path.abspath(newcase)
        newcasename = os.path.basename(newcase)
        os.makedirs(newcaseroot)
        clone = self.copy(newcasename = newcasename, newcaseroot = newcaseroot)

        return clone

    def flush(self):
        pass

# ========================================================================
# Tests of CaseFake
# ========================================================================

class TestCaseFake(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_create_clone(self):
        # Setup
        old_caseroot = os.path.join(self.tempdir, 'oldcase')
        oldcase = CaseFake(old_caseroot)
        oldcase.set_value('foo', 'bar')

        # Exercise
        new_caseroot = os.path.join(self.tempdir, 'newcase')
        clone = oldcase.create_clone(new_caseroot)

        # Verify
        self.assertEqual('bar', clone.get_value('foo'))
        self.assertEqual('newcase', clone.get_value('CASE'))
        self.assertEqual('newcase', clone.get_value('CASEBASEID'))
        self.assertEqual(new_caseroot, clone.get_value('CASEROOT'))

# ========================================================================
# Fake version of SystemTestsCompareTwo that overrides some functionality for
# the sake of unit testing
# ========================================================================

class SystemTestsCompareTwoFake(SystemTestsCompareTwoClone):
    def __init__(self,
                 case1,
                 run_two_suffix = 'test',
                 case2setup_raises_exception=False):
        """
        Initialize a SystemTestsCompareTwoFake object

        Args:
            case1 (CaseFake): existing case
            run_two_suffix (str, optional): Suffix used for the second run. Defaults to 'test'.
            case2setup_raises_exception (bool, optional): If True, then the call
                to _case_two_setup will raise an exception. Defaults to False.
        """

        self._case2setup_raises_exception = case2setup_raises_exception

        SystemTestsCompareTwoClone.__init__(
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
    # Fake implementations of methods that are typically provided by
    # SystemTestsCompareTwoClone
    #
    # Since we're overriding these, their functionality is untested!
    # ------------------------------------------------------------------------

    def _case_from_existing_caseroot(self, caseroot):
        """
        Returns a CaseFake object instead of a Case object
        """
        return CaseFake(caseroot, create_case_root=False)

    # ------------------------------------------------------------------------
    # Fake implementations of methods that are typically provided by the
    # individual test
    #
    # The values set here are asserted against in some unit tests
    # ------------------------------------------------------------------------

    def _common_setup(self):
        self._case.set_value('var_set_in_common_setup', 'common_val')

    def _case_one_setup(self):
        self._case.set_value('var_set_in_setup', 'case1val')

    def _case_two_setup(self):
        self._case.set_value('var_set_in_setup', 'case2val')
        if self._case2setup_raises_exception:
            raise RuntimeError

# ========================================================================
# Test class itself
# ========================================================================

class TestSystemTestsCompareTwoClone(unittest.TestCase):

    def setUp(self):
        # create a sandbox in which case directories can be created
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_setup(self):
        # Ensure that test setup properly sets up case 1 and case 2

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)

        # Exercise
        mytest = SystemTestsCompareTwoFake(case1)

        # Verify
        self.assertEqual('common_val',
                         mytest._case1.get_value('var_set_in_common_setup'))
        self.assertEqual('common_val',
                         mytest._case2.get_value('var_set_in_common_setup'))
        self.assertEqual('case1val',
                         mytest._case1.get_value('var_set_in_setup'))
        self.assertEqual('case2val',
                         mytest._case2.get_value('var_set_in_setup'))

    def test_setup_case2_exists(self):
        # If case2 already exists, then setup code should not be called

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        os.makedirs(os.path.join(case1root, 'case1.test'))

        # Exercise
        mytest = SystemTestsCompareTwoFake(case1,
                                           run_two_suffix = 'test')

        # Verify: variables set in various setup methods should not be set
        # (In the real world - i.e., outside of this unit testing fakery - these
        # values would be set when the Case objects are created.)
        self.assertIsNone(mytest._case1.get_value('var_set_in_common_setup'))
        self.assertIsNone(mytest._case2.get_value('var_set_in_common_setup'))
        self.assertIsNone(mytest._case1.get_value('var_set_in_setup'))
        self.assertIsNone(mytest._case2.get_value('var_set_in_setup'))

    def test_setup_error(self):
        # If there is an error in setup, an exception should be raised and the
        # case2 directory should be removed

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)

        # Exercise
        with self.assertRaises(Exception):
            mytest = SystemTestsCompareTwoFake(case1,
                                               run_two_suffix = 'test',
                                               case2setup_raises_exception = True)

        # Verify
        self.assertFalse(os.path.exists(os.path.join(case1root, 'case1.test')))
