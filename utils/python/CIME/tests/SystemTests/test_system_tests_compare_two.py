#!/usr/bin/env python

"""
This module contains unit tests of the core logic in SystemTestsCompareTwo.
"""

# Ignore privacy concerns for unit tests, so that unit tests can access
# protected members of the system under test
#
# pylint:disable=protected-access

import unittest
from collections import namedtuple
import os
import shutil
import tempfile

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
import CIME.test_status as test_status
from CIME.tests.case_fake import CaseFake

# ========================================================================
# Structure for storing information about calls made to methods
# ========================================================================

# You can create a Call object to record a single call made to a method:
#
# Call(method, arguments)
#     method (str): name of method
#     arguments (dict): dictionary mapping argument names to values
#
# Example:
#     If you want to record a call to foo(bar = 1, baz = 2):
#         somecall = Call(method = 'foo', arguments = {'bar': 1, 'baz': 2})
#     Or simply:
#         somecall = Call('foo', {'bar': 1, 'baz': 2})
Call = namedtuple('Call', ['method', 'arguments'])

# ========================================================================
# Names of methods for which we want to record calls
# ========================================================================

# We use constants for these method names because, in some cases, a typo in a
# hard-coded string could cause a test to always pass, which would be a Bad
# Thing.
#
# For now the names of the constants match the strings they equate to, which
# match the actual method names. But it's fine if this doesn't remain the case
# moving forward (which is another reason to use constants rather than
# hard-coded strings in the tests).

METHOD_component_compare_test = "_component_compare_test"
METHOD_link_to_case2_output = "_link_to_case2_output"
METHOD_run_indv = "_run_indv"

# ========================================================================
# Fake version of SystemTestsCompareTwo that overrides some functionality for
# the sake of unit testing
# ========================================================================

# A SystemTestsCompareTwoFake object can be controlled to fail at a given
# point. See the documentation in its __init__ method for details.
#
# It logs what stubbed-out methods have been called in its log attribute; this
# is a list of Call objects (see above for their definition).

class SystemTestsCompareTwoFake(SystemTestsCompareTwo):
    def __init__(self,
                 case1,
                 run_one_suffix = 'base',
                 run_two_suffix = 'test',
                 case2setup_raises_exception = False,
                 run_one_should_pass = True,
                 run_two_should_pass = True):
        """
        Initialize a SystemTestsCompareTwoFake object

        Args:
            case1 (CaseFake): existing case
            run_one_suffix (str, optional): Suffix used for first run. Defaults
                to 'base'. Currently MUST be 'base'.
            run_two_suffix (str, optional): Suffix used for the second run. Defaults to 'test'.
            case2setup_raises_exception (bool, optional): If True, then the call
                to _case_two_setup will raise an exception. Default is False.
            run_one_should_pass (bool, optional): Whether the run_indv method should
                pass for the first run. Default is True, meaning it will pass.
            run_two_should_pass (bool, optional): Whether the run_indv method should
                pass for the second run. Default is True, meaning it will pass.
        """

        self._case2setup_raises_exception = case2setup_raises_exception

        # NOTE(wjs, 2016-08-03) Currently, due to limitations in the test
        # infrastructure, run_one_suffix MUST be 'base'. However, I'm keeping it
        # as an explicit argument to the constructor so that it's easy to relax
        # this requirement later: To relax this assumption, remove the following
        # assertion and add run_one_suffix as an argument to
        # SystemTestsCompareTwo.__init__
        assert(run_one_suffix == 'base')

        SystemTestsCompareTwo.__init__(
            self,
            case1,
            separate_builds = False,
            run_two_suffix = run_two_suffix)

        # Need to tell test status that case has been built, since this is
        # checked in the run call
        with self._test_status:
            for phase in test_status.CORE_PHASES[:-1]:
                self._test_status.set_status(phase, test_status.TEST_PASS_STATUS)

        self.run_pass_casenames = []
        if run_one_should_pass:
            self.run_pass_casenames.append(self._case1.get_value('CASE'))
        if run_two_should_pass:
            self.run_pass_casenames.append(self._case2.get_value('CASE'))

        self.log = []

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
    # SystemTestsCommon
    # ------------------------------------------------------------------------

    def run_indv(self, suffix="base", st_archive=False):
        """
        This fake implementation appends to the log and raises an exception if
        it's supposed to

        Note that the Call object appended to the log has the current CASE name
        in addition to the method arguments. (This is mainly to ensure that the
        proper suffix is used for the proper case, but this extra check can be
        removed if it's a maintenance problem.)
        """
        casename = self._case.get_value('CASE')
        self.log.append(Call(METHOD_run_indv,
                             {'suffix': suffix, 'CASE': casename}))

        # Determine whether we should raise an exception
        #
        # It's important that this check be based on some attribute of the
        # self._case object, to ensure that the right case has been activated
        # for this call to run_indv (e.g., to catch if we forgot to activate
        # case2 before the second call to run_indv).
        if casename not in self.run_pass_casenames:
            raise RuntimeError('casename not in run_pass_casenames')

    def _component_compare_test(self, suffix1, suffix2):
        # Trying to use the real version of _component_compare_test would pull
        # too much baggage into these tests. Since the return value from this
        # method isn't important, it's sufficient for the tests of this class to
        # just ensure that _component_compare_test was actually called
        # correctly.
        #
        # An alternative would be to extract the main work of
        # _component_compare_test into a different method that returns a True
        # (success) / False (failure) result, with _component_compare_test then
        # updating test_status appropriately. Then we could override that new
        # method in this Fake class, using the true implementation of
        # _component_compare_test. Then the test verification would include
        # verification that TestStatus is set correctly for the COMPARE
        # phase. But that seems more about testing _component_compare_test than
        # testing SystemTestsCompareTwo itself, so I don't see much added value
        # of that.

        self.log.append(Call(METHOD_component_compare_test,
                             {'suffix1': suffix1, 'suffix2': suffix2}))

    def _check_for_memleak(self):
        pass

    # ------------------------------------------------------------------------
    # Fake implementations of methods that are typically provided by
    # SystemTestsCompareTwo
    #
    # Since we're overriding these, their functionality is untested here!
    # (Though note that _link_to_case2_output is tested elsewhere.)
    # ------------------------------------------------------------------------

    def _case_from_existing_caseroot(self, caseroot):
        """
        Returns a CaseFake object instead of a Case object
        """
        return CaseFake(caseroot, create_case_root=False)

    def _link_to_case2_output(self):
        self.log.append(Call(METHOD_link_to_case2_output, {}))

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

class TestSystemTestsCompareTwo(unittest.TestCase):

    def setUp(self):
        self.original_wd = os.getcwd()
        # create a sandbox in which case directories can be created
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        # Some tests trigger a chdir call in the SUT; make sure we return to the
        # original directory at the end of the test
        os.chdir(self.original_wd)

        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_setup(self):
        # Ensure that test setup properly sets up case 1 and case 2

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        case1.set_value('var_preset', 'preset_value')

        # Exercise
        mytest = SystemTestsCompareTwoFake(case1)

        # Verify
        # Make sure that pre-existing values in case1 are copied to case2 (via
        # clone)
        self.assertEqual('preset_value',
                         mytest._case2.get_value('var_preset'))

        # Make sure that _common_setup is called for both
        self.assertEqual('common_val',
                         mytest._case1.get_value('var_set_in_common_setup'))
        self.assertEqual('common_val',
                         mytest._case2.get_value('var_set_in_common_setup'))

        # Make sure that _case_one_setup and _case_two_setup are called
        # appropriately
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

        # Verify:

        # Make sure that case2 object is set (i.e., that it doesn't remain None)
        self.assertEqual('case1.test', mytest._case2.get_value('CASE'))

        # Variables set in various setup methods should not be set
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
            SystemTestsCompareTwoFake(case1,
                                      run_two_suffix = 'test',
                                      case2setup_raises_exception = True)

        # Verify
        self.assertFalse(os.path.exists(os.path.join(case1root, 'case1.test')))

    def test_run_phase_passes(self):
        # Make sure the run phase behaves properly when all runs succeed.

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1)

        # Exercise
        mytest.run()

        # Verify
        self.assertEqual(test_status.TEST_PASS_STATUS,
                         mytest._test_status.get_status(test_status.RUN_PHASE))

    def test_run_phase_internal_calls(self):
        # Make sure that the correct calls are made to methods stubbed out by
        # SystemTestsCompareTwoFake (when runs succeed)
        #
        # The point of this is: A number of methods called from the run_phase
        # method are stubbed out in the Fake test implementation, because their
        # actions are awkward in these unit tests. But we still want to make
        # sure that those methods actually got called correctly.

        # Setup
        run_one_suffix = 'base'
        run_two_suffix = 'run2'
        casename = 'mytest'
        case1root = os.path.join(self.tempdir, casename)
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1,
            run_one_suffix = run_one_suffix,
            run_two_suffix = run_two_suffix)

        # Exercise
        mytest.run()

        # Verify
        expected_calls = [
            Call(METHOD_run_indv,
                 {'suffix': run_one_suffix, 'CASE': casename}),
            Call(METHOD_run_indv,
                 {'suffix': run_two_suffix, 'CASE': '%s.%s'%(casename, run_two_suffix)}),
            Call(METHOD_link_to_case2_output, {}),
            Call(METHOD_component_compare_test,
                 {'suffix1': run_one_suffix, 'suffix2': run_two_suffix})]
        self.assertEqual(expected_calls, mytest.log)

    def test_run1_fails(self):
        # Make sure that a failure in run1 is reported correctly

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1,
                                           run_one_should_pass = False)

        # Exercise
        mytest.run()

        # Verify
        self.assertEqual(test_status.TEST_FAIL_STATUS,
                         mytest._test_status.get_status(test_status.RUN_PHASE))

    def test_run2_fails(self):
        # Make sure that a failure in run2 is reported correctly

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1,
                                           run_two_should_pass = False)

        # Exercise
        mytest.run()

        # Verify
        self.assertEqual(test_status.TEST_FAIL_STATUS,
                         mytest._test_status.get_status(test_status.RUN_PHASE))

