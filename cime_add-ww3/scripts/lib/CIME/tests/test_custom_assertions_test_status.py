#!/usr/bin/env python

"""
This module contains unit tests of CustomAssertionsTestStatus
"""

import unittest
from CIME import test_status
from CIME.tests.custom_assertions_test_status import CustomAssertionsTestStatus

class TestCustomAssertions(CustomAssertionsTestStatus):

    _UNEXPECTED_COMMENT = test_status.TEST_UNEXPECTED_FAILURE_COMMENT_START + ' blah)'

    @staticmethod
    def output_line(status, test_name, phase, extra=''):
        output = status + ' ' + test_name + ' ' + phase
        if extra:
            output += ' ' + extra
        output += '\n'
        return output

    def test_assertPhaseAbsent_passes(self):
        """assert_phase_absent should pass when the phase is absent for
        the given test_name"""
        test_name1 = 'my.test.name1'
        test_name2 = 'my.test.name2'
        output = self.output_line('PASS', test_name1, 'PHASE1')
        output += self.output_line('PASS', test_name2, 'PHASE2')

        self.assert_phase_absent(output, 'PHASE2', test_name1)
        self.assert_phase_absent(output, 'PHASE1', test_name2)

    def test_assertPhaseAbsent_fails(self):
        """assert_phase_absent should fail when the phase is present for
        the given test_name"""
        test_name = 'my.test.name'
        output = self.output_line('PASS', test_name, 'PHASE1')

        with self.assertRaises(AssertionError):
            self.assert_phase_absent(output, 'PHASE1', test_name)

    def test_assertCorePhases_passes(self):
        """assert_core_phases passes when it should"""
        output = ''
        fails = [test_status.CORE_PHASES[1]]
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            if phase in fails:
                status = test_status.TEST_FAIL_STATUS
            else:
                status = test_status.TEST_PASS_STATUS
            output = output + self.output_line(status, test_name, phase)

        self.assert_core_phases(output, test_name, fails)

    def test_assertCorePhases_missingPhase_fails(self):
        """assert_core_phases fails if there is a missing phase"""
        output = ''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            if phase != test_status.CORE_PHASES[1]:
                output = output + self.output_line(test_status.TEST_PASS_STATUS, test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, test_name, fails=[])

    def test_assertCorePhases_wrongStatus_fails(self):
        """assert_core_phases fails if a phase has the wrong status"""
        output =''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            output = output + self.output_line(test_status.TEST_PASS_STATUS, test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, test_name, fails=[test_status.CORE_PHASES[1]])

    def test_assertCorePhases_wrongName_fails(self):
        """assert_core_phases fails if the test name is wrong"""
        output =''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            output = output + self.output_line(test_status.TEST_PASS_STATUS, test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, 'my.test', fails=[])

    # Note: Basic functionality of assert_status_of_phase is covered sufficiently via
    # tests of assert_core_phases. Below we just cover some other aspects that aren't
    # already covered.

    def test_assertStatusOfPhase_withExtra_passes(self):
        """Make sure assert_status_of_phase passes when there is some extra text at the
        end of the line"""
        test_name = 'my.test.name'
        output = self.output_line(test_status.TEST_FAIL_STATUS,
                                  test_name,
                                  test_status.CORE_PHASES[0],
                                  extra=test_status.TEST_EXPECTED_FAILURE_COMMENT)
        self.assert_status_of_phase(output,
                                    test_status.TEST_FAIL_STATUS,
                                    test_status.CORE_PHASES[0],
                                    test_name)

    def test_assertStatusOfPhase_xfailNo_passes(self):
        """assert_status_of_phase should pass when xfail='no' and there is no
        EXPECTED/UNEXPECTED on the line"""
        test_name = 'my.test.name'
        output = self.output_line(test_status.TEST_FAIL_STATUS,
                                  test_name,
                                  test_status.CORE_PHASES[0])
        self.assert_status_of_phase(output,
                                    test_status.TEST_FAIL_STATUS,
                                    test_status.CORE_PHASES[0],
                                    test_name,
                                    xfail='no')
        # While we're at it, also test assert_num_expected_unexpected_fails
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=0,
                                                  num_unexpected=0)

    def test_assertStatusOfPhase_xfailNo_fails(self):
        """assert_status_of_phase should fail when xfail='no' but the line contains the
        EXPECTED comment"""
        test_name = 'my.test.name'
        output = self.output_line(test_status.TEST_FAIL_STATUS,
                                  test_name,
                                  test_status.CORE_PHASES[0],
                                  extra=test_status.TEST_EXPECTED_FAILURE_COMMENT)

        with self.assertRaises(AssertionError):
            self.assert_status_of_phase(output,
                                        test_status.TEST_FAIL_STATUS,
                                        test_status.CORE_PHASES[0],
                                        test_name,
                                        xfail='no')
        # While we're at it, also test assert_num_expected_unexpected_fails
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=1,
                                                  num_unexpected=0)

    def test_assertStatusOfPhase_xfailExpected_passes(self):
        """assert_status_of_phase should pass when xfail='expected' and the line contains
        the EXPECTED comment"""
        test_name = 'my.test.name'
        output = self.output_line(test_status.TEST_FAIL_STATUS,
                                  test_name,
                                  test_status.CORE_PHASES[0],
                                  extra=test_status.TEST_EXPECTED_FAILURE_COMMENT)
        self.assert_status_of_phase(output,
                                    test_status.TEST_FAIL_STATUS,
                                    test_status.CORE_PHASES[0],
                                    test_name,
                                    xfail='expected')
        # While we're at it, also test assert_num_expected_unexpected_fails
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=1,
                                                  num_unexpected=0)

    def test_assertStatusOfPhase_xfailExpected_fails(self):
        """assert_status_of_phase should fail when xfail='expected' but the line does NOT contain
        the EXPECTED comment"""
        test_name = 'my.test.name'
        # Note that the line contains the UNEXPECTED comment, but not the EXPECTED comment
        # (we assume that if the assertion correctly fails in this case, then it will also
        # correctly handle the case where neither the EXPECTED nor UNEXPECTED comment is
        # present).
        output = self.output_line(test_status.TEST_FAIL_STATUS,
                                  test_name,
                                  test_status.CORE_PHASES[0],
                                  extra=self._UNEXPECTED_COMMENT)

        with self.assertRaises(AssertionError):
            self.assert_status_of_phase(output,
                                        test_status.TEST_FAIL_STATUS,
                                        test_status.CORE_PHASES[0],
                                        test_name,
                                        xfail='expected')
        # While we're at it, also test assert_num_expected_unexpected_fails
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=0,
                                                  num_unexpected=1)

    def test_assertStatusOfPhase_xfailUnexpected_passes(self):
        """assert_status_of_phase should pass when xfail='unexpected' and the line contains
        the UNEXPECTED comment"""
        test_name = 'my.test.name'
        output = self.output_line(test_status.TEST_FAIL_STATUS,
                                  test_name,
                                  test_status.CORE_PHASES[0],
                                  extra=self._UNEXPECTED_COMMENT)
        self.assert_status_of_phase(output,
                                    test_status.TEST_FAIL_STATUS,
                                    test_status.CORE_PHASES[0],
                                    test_name,
                                    xfail='unexpected')
        # While we're at it, also test assert_num_expected_unexpected_fails
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=0,
                                                  num_unexpected=1)

    def test_assertStatusOfPhase_xfailUnexpected_fails(self):
        """assert_status_of_phase should fail when xfail='unexpected' but the line does NOT
        contain the UNEXPECTED comment"""
        test_name = 'my.test.name'
        # Note that the line contains the EXPECTED comment, but not the UNEXPECTED comment
        # (we assume that if the assertion correctly fails in this case, then it will also
        # correctly handle the case where neither the EXPECTED nor UNEXPECTED comment is
        # present).
        output = self.output_line(test_status.TEST_FAIL_STATUS,
                                  test_name,
                                  test_status.CORE_PHASES[0],
                                  extra=test_status.TEST_EXPECTED_FAILURE_COMMENT)

        with self.assertRaises(AssertionError):
            self.assert_status_of_phase(output,
                                        test_status.TEST_FAIL_STATUS,
                                        test_status.CORE_PHASES[0],
                                        test_name,
                                        xfail='unexpected')
        # While we're at it, also test assert_num_expected_unexpected_fails
        self.assert_num_expected_unexpected_fails(output,
                                                  num_expected=1,
                                                  num_unexpected=0)

if __name__ == '__main__':
    unittest.main()
