#!/usr/bin/env python

import unittest
import os
from CIME import test_status
from CIME import expected_fails
from CIME.tests.custom_assertions_test_status import CustomAssertionsTestStatus

class TestTestStatus(CustomAssertionsTestStatus):

    _TESTNAME = 'fake_test'

    # An arbitrary phase we can use when we want to work with a non-core phase
    _NON_CORE_PHASE = test_status.MEMLEAK_PHASE

    def setUp(self):
        self._ts = test_status.TestStatus(test_dir=os.path.join('nonexistent', 'path'),
                                          test_name=self._TESTNAME,
                                          no_io=True)
        self._set_core_phases_to_pass()

    def _set_core_phases_to_pass(self):
        """Set all core phases of self._ts to pass status"""
        with self._ts:
            for phase in test_status.CORE_PHASES:
                self._ts.set_status(phase, test_status.TEST_PASS_STATUS)

    def _set_last_core_phase_to_fail(self):
        """Sets the last core phase to FAIL

        Returns the name of this phase"""
        fail_phase = test_status.CORE_PHASES[-1]
        self._set_phase_to_status(fail_phase, test_status.TEST_FAIL_STATUS)
        return fail_phase

    def _set_phase_to_status(self, phase, status):
        """Set given phase to given status"""
        with self._ts:
            self._ts.set_status(phase, status)

    # ------------------------------------------------------------------------
    # Tests of TestStatus.phase_statuses_dump
    # ------------------------------------------------------------------------

    def test_psdump_corePhasesPass(self):
        output = self._ts.phase_statuses_dump()
        self.assert_core_phases(output, self._TESTNAME, fails=[])
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=0,
                                                  num_unexpected=0)

    def test_psdump_oneCorePhaseFails(self):
        fail_phase = self._set_last_core_phase_to_fail()
        output = self._ts.phase_statuses_dump()
        self.assert_core_phases(output, self._TESTNAME, fails=[fail_phase])
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=0,
                                                  num_unexpected=0)

    def test_psdump_oneCorePhaseFailsAbsentFromXFails(self):
        """One phase fails. There is an expected fails list, but that phase is not in it."""
        fail_phase = self._set_last_core_phase_to_fail()
        xfails = expected_fails.ExpectedFails()
        xfails.add_failure(phase=self._NON_CORE_PHASE,
                           expected_status=test_status.TEST_FAIL_STATUS)
        output = self._ts.phase_statuses_dump(xfails=xfails)
        self.assert_status_of_phase(output,
                                    test_status.TEST_FAIL_STATUS,
                                    fail_phase,
                                    self._TESTNAME,
                                    xfail='no')
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=0,
                                                  num_unexpected=0)

    def test_psdump_oneCorePhaseFailsInXFails(self):
        """One phase fails. That phase is in the expected fails list."""
        fail_phase = self._set_last_core_phase_to_fail()
        xfails = expected_fails.ExpectedFails()
        xfails.add_failure(phase=fail_phase,
                           expected_status=test_status.TEST_FAIL_STATUS)
        output = self._ts.phase_statuses_dump(xfails=xfails)
        self.assert_status_of_phase(output,
                                    test_status.TEST_FAIL_STATUS,
                                    fail_phase,
                                    self._TESTNAME,
                                    xfail='expected')
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=1,
                                                  num_unexpected=0)

    def test_psdump_oneCorePhasePassesInXFails(self):
        """One phase passes despite being in the expected fails list."""
        xfail_phase = test_status.CORE_PHASES[-1]
        xfails = expected_fails.ExpectedFails()
        xfails.add_failure(phase=xfail_phase,
                           expected_status=test_status.TEST_FAIL_STATUS)
        output = self._ts.phase_statuses_dump(xfails=xfails)
        self.assert_status_of_phase(output,
                                    test_status.TEST_PASS_STATUS,
                                    xfail_phase,
                                    self._TESTNAME,
                                    xfail='unexpected')
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=0,
                                                  num_unexpected=1)

    def test_psdump_skipPasses(self):
        """With the skip_passes argument, only non-passes should appear"""
        fail_phase = self._set_last_core_phase_to_fail()
        output = self._ts.phase_statuses_dump(skip_passes=True)
        self.assert_status_of_phase(output,
                                    test_status.TEST_FAIL_STATUS,
                                    fail_phase,
                                    self._TESTNAME,
                                    xfail='no')
        for phase in test_status.CORE_PHASES:
            if phase != fail_phase:
                self.assert_phase_absent(output, phase, self._TESTNAME)

    def test_psdump_unexpectedPass_shouldBePresent(self):
        """Even with the skip_passes argument, an unexpected PASS should be present"""
        xfail_phase = test_status.CORE_PHASES[-1]
        xfails = expected_fails.ExpectedFails()
        xfails.add_failure(phase=xfail_phase,
                           expected_status=test_status.TEST_FAIL_STATUS)
        output = self._ts.phase_statuses_dump(skip_passes=True, xfails=xfails)
        self.assert_status_of_phase(output,
                                    test_status.TEST_PASS_STATUS,
                                    xfail_phase,
                                    self._TESTNAME,
                                    xfail='unexpected')
        for phase in test_status.CORE_PHASES:
            if phase != xfail_phase:
                self.assert_phase_absent(output, phase, self._TESTNAME)

if __name__ == '__main__':
    unittest.main()
