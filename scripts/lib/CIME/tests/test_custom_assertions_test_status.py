#!/usr/bin/env python

"""
This module contains unit tests of CustomAssertionsTestStatus
"""

import unittest
from CIME import test_status
from CIME.tests.custom_assertions_test_status import CustomAssertionsTestStatus

class TestCustomAssertions(CustomAssertionsTestStatus):

    @staticmethod
    def output_line(status, test_name, phase):
        return status + ' ' + test_name + ' ' + phase + '\n'

    # Note: No explicit tests of assert_status_of_phase: this is covered
    # sufficiently via tests of assert_core_phases

    def test_assert_phase_absent_passes(self):
        """assert_phase_absent should pass when the phase is absent for
        the given test_name"""
        test_name1 = 'my.test.name1'
        test_name2 = 'my.test.name2'
        output = self.output_line('PASS', test_name1, 'PHASE1')
        output += self.output_line('PASS', test_name2, 'PHASE2')

        self.assert_phase_absent(output, 'PHASE2', test_name1)
        self.assert_phase_absent(output, 'PHASE1', test_name2)

    def test_assert_phase_absent_fails(self):
        """assert_phase_absent should fail when the phase is present for
        the given test_name"""
        test_name = 'my.test.name'
        output = self.output_line('PASS', test_name, 'PHASE1')

        with self.assertRaises(AssertionError):
            self.assert_phase_absent(output, 'PHASE1', test_name)

    def test_assert_core_phases_passes(self):
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

    def test_assert_core_phases_missing_phase_fails(self):
        """assert_core_phases fails if there is a missing phase"""
        output = ''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            if phase != test_status.CORE_PHASES[1]:
                output = output + self.output_line(test_status.TEST_PASS_STATUS, test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, test_name, fails=[])

    def test_assert_core_phases_wrong_status_fails(self):
        """assert_core_phases fails if a phase has the wrong status"""
        output =''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            output = output + self.output_line(test_status.TEST_PASS_STATUS, test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, test_name, fails=[test_status.CORE_PHASES[1]])

    def test_assert_core_phases_wrong_name_fails(self):
        """assert_core_phases fails if the test name is wrong"""
        output =''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            output = output + self.output_line(test_status.TEST_PASS_STATUS, test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, 'my.test', fails=[])

if __name__ == '__main__':
    unittest.main()
