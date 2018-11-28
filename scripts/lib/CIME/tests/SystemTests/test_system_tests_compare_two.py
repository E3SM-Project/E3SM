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

METHOD_case_one_custom_prerun_action = "_case_one_custom_prerun_action"
METHOD_case_one_custom_postrun_action = "_case_one_custom_postrun_action"
METHOD_case_two_custom_prerun_action = "_case_two_custom_prerun_action"
METHOD_case_two_custom_postrun_action = "_case_two_custom_postrun_action"
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
                 separate_builds = False,
                 multisubmit = False,
                 case2setup_raises_exception = False,
                 run_one_should_pass = True,
                 run_two_should_pass = True,
                 compare_should_pass = True):
        """
        Initialize a SystemTestsCompareTwoFake object

        The core test phases prior to RUN_PHASE are set to TEST_PASS_STATUS;
        RUN_PHASE is left unset (as is any later phase)

        Args:
            case1 (CaseFake): existing case
            run_one_suffix (str, optional): Suffix used for first run. Defaults
                to 'base'. Currently MUST be 'base'.
            run_two_suffix (str, optional): Suffix used for the second run. Defaults to 'test'.
            separate_builds (bool, optional): Passed to SystemTestsCompareTwo.__init__
            multisubmit (bool, optional): Passed to SystemTestsCompareTwo.__init__
            case2setup_raises_exception (bool, optional): If True, then the call
                to _case_two_setup will raise an exception. Default is False.
            run_one_should_pass (bool, optional): Whether the run_indv method should
                pass for the first run. Default is True, meaning it will pass.
            run_two_should_pass (bool, optional): Whether the run_indv method should
                pass for the second run. Default is True, meaning it will pass.
            compare_should_pass (bool, optional): Whether the comparison between the two
                cases should pass. Default is True, meaning it will pass.
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
            separate_builds = separate_builds,
            run_two_suffix = run_two_suffix,
            multisubmit = multisubmit)

        # Need to tell test status that all phases prior to the run phase have
        # passed, since this is checked in the run call (at least for the build
        # phase status)
        with self._test_status:
            for phase in test_status.CORE_PHASES:
                if phase == test_status.RUN_PHASE:
                    break
                self._test_status.set_status(phase, test_status.TEST_PASS_STATUS)

        self.run_pass_caseroot = []
        if run_one_should_pass:
            self.run_pass_caseroot.append(self._case1.get_value('CASEROOT'))
        if run_two_should_pass:
            self.run_pass_caseroot.append(self._case2.get_value('CASEROOT'))

        self.compare_should_pass = compare_should_pass

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
        caseroot = self._case.get_value('CASEROOT')
        self.log.append(Call(METHOD_run_indv,
                             {'suffix': suffix, 'CASEROOT': caseroot}))

        # Determine whether we should raise an exception
        #
        # It's important that this check be based on some attribute of the
        # self._case object, to ensure that the right case has been activated
        # for this call to run_indv (e.g., to catch if we forgot to activate
        # case2 before the second call to run_indv).
        if caseroot not in self.run_pass_caseroot:
            raise RuntimeError('caseroot not in run_pass_caseroot')

    def _do_compare_test(self, suffix1, suffix2):
        """
        This fake implementation allows controlling whether compare_test
        passes or fails
        """
        return (self.compare_should_pass, "no comment")

    def _check_for_memleak(self):
        pass

    def _st_archive_case_test(self):
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

    def _case_one_custom_prerun_action(self):
        self.log.append(Call(METHOD_case_one_custom_prerun_action, {}))

    def _case_one_custom_postrun_action(self):
        self.log.append(Call(METHOD_case_one_custom_postrun_action, {}))

    def _case_two_custom_prerun_action(self):
        self.log.append(Call(METHOD_case_two_custom_prerun_action, {}))

    def _case_two_custom_postrun_action(self):
        self.log.append(Call(METHOD_case_two_custom_postrun_action, {}))

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

    def get_caseroots(self, casename='mytest'):
        """
        Returns a tuple (case1root, case2root)
        """
        case1root = os.path.join(self.tempdir, casename)
        case2root = os.path.join(case1root, 'case2', casename)
        return case1root, case2root

    def get_compare_phase_name(self, mytest):
        """
        Returns a string giving the compare phase name for this test
        """
        run_one_suffix = mytest._run_one_suffix
        run_two_suffix = mytest._run_two_suffix
        compare_phase_name = "{}_{}_{}".format(test_status.COMPARE_PHASE,
                                               run_one_suffix,
                                               run_two_suffix)
        return compare_phase_name

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

    def test_setup_separate_builds_sharedlibroot(self):
        # If we're using separate_builds, the two cases should still use
        # the same sharedlibroot

        # Setup
        case1root, _ = self.get_caseroots()
        case1 = CaseFake(case1root)
        case1.set_value("SHAREDLIBROOT", os.path.join(case1root, "sharedlibroot"))

        # Exercise
        mytest = SystemTestsCompareTwoFake(case1,
                                           separate_builds = True)

        # Verify
        self.assertEqual(case1.get_value("SHAREDLIBROOT"),
                         mytest._case2.get_value("SHAREDLIBROOT"))

    def test_setup_case2_exists(self):
        # If case2 already exists, then setup code should not be called

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        os.makedirs(os.path.join(case1root, 'case2','case1'))

        # Exercise
        mytest = SystemTestsCompareTwoFake(case1,
                                           run_two_suffix = 'test')

        # Verify:

        # Make sure that case2 object is set (i.e., that it doesn't remain None)
        self.assertEqual('case1', mytest._case2.get_value('CASE'))

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
        case1root, case2root = self.get_caseroots()
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1,
            run_one_suffix = run_one_suffix,
            run_two_suffix = run_two_suffix)

        # Exercise
        mytest.run()

        # Verify
        expected_calls = [
            Call(METHOD_case_one_custom_prerun_action, {}),
            Call(METHOD_run_indv,
                 {'suffix': run_one_suffix, 'CASEROOT': case1root}),
            Call(METHOD_case_one_custom_postrun_action, {}),
            Call(METHOD_case_two_custom_prerun_action, {}),
            Call(METHOD_run_indv,
                 {'suffix': run_two_suffix, 'CASEROOT': case2root}),
            Call(METHOD_case_two_custom_postrun_action, {}),
            Call(METHOD_link_to_case2_output, {})
        ]
        self.assertEqual(expected_calls, mytest.log)

    def test_run_phase_internal_calls_multisubmit_phase1(self):
        # Make sure that the correct calls are made to methods stubbed out by
        # SystemTestsCompareTwoFake (when runs succeed), when we have a
        # multi-submit test, in the first phase

        # Setup
        run_one_suffix = 'base'
        run_two_suffix = 'run2'
        case1root, _ = self.get_caseroots()
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(
            case1 = case1,
            run_one_suffix = run_one_suffix,
            run_two_suffix = run_two_suffix,
            multisubmit = True)
        # RESUBMIT=1 signals first phase
        case1.set_value("RESUBMIT", 1)

        # Exercise
        mytest.run()

        # Verify
        expected_calls = [
            Call(METHOD_case_one_custom_prerun_action, {}),
            Call(METHOD_run_indv,
                 {'suffix': run_one_suffix, 'CASEROOT': case1root}),
            Call(METHOD_case_one_custom_postrun_action, {}),
        ]
        self.assertEqual(expected_calls, mytest.log)

        # Also verify that comparison is NOT called:
        compare_phase_name = self.get_compare_phase_name(mytest)
        self.assertEqual(test_status.TEST_PEND_STATUS, mytest._test_status.get_status(compare_phase_name))

    def test_run_phase_internal_calls_multisubmit_phase2(self):
        # Make sure that the correct calls are made to methods stubbed out by
        # SystemTestsCompareTwoFake (when runs succeed), when we have a
        # multi-submit test, in the second phase

        # Setup
        run_one_suffix = 'base'
        run_two_suffix = 'run2'
        case1root, case2root = self.get_caseroots()
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(
            case1 = case1,
            run_one_suffix = run_one_suffix,
            run_two_suffix = run_two_suffix,
            multisubmit = True,
            compare_should_pass = True)
        # RESUBMIT=0 signals second phase
        case1.set_value("RESUBMIT", 0)

        # Exercise
        mytest.run()

        # Verify
        expected_calls = [
            Call(METHOD_case_two_custom_prerun_action, {}),
            Call(METHOD_run_indv,
                 {'suffix': run_two_suffix, 'CASEROOT': case2root}),
            Call(METHOD_case_two_custom_postrun_action, {}),
            Call(METHOD_link_to_case2_output, {})
        ]
        self.assertEqual(expected_calls, mytest.log)

        # Also verify that comparison is called:
        compare_phase_name = self.get_compare_phase_name(mytest)
        self.assertEqual(test_status.TEST_PASS_STATUS,
                         mytest._test_status.get_status(compare_phase_name))

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

    def test_compare_passes(self):
        # Make sure that a pass in the comparison is reported correctly

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1,
                                           compare_should_pass = True)

        # Exercise
        mytest.run()

        # Verify
        compare_phase_name = self.get_compare_phase_name(mytest)
        self.assertEqual(test_status.TEST_PASS_STATUS,
                         mytest._test_status.get_status(compare_phase_name))

    def test_compare_fails(self):
        # Make sure that a failure in the comparison is reported correctly

        # Setup
        case1root = os.path.join(self.tempdir, 'case1')
        case1 = CaseFake(case1root)
        mytest = SystemTestsCompareTwoFake(case1,
                                           compare_should_pass = False)

        # Exercise
        mytest.run()

        # Verify
        compare_phase_name = self.get_compare_phase_name(mytest)
        self.assertEqual(test_status.TEST_FAIL_STATUS,
                         mytest._test_status.get_status(compare_phase_name))

if __name__ == "__main__":
    unittest.main(verbosity=2, catchbreak=True)
