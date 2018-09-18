#!/usr/bin/env python

import unittest
import shutil
import os
import tempfile
import re
import six
import six_additions
from CIME.cs_status import cs_status
from CIME import test_status

# ========================================================================
# Custom assertions
# ========================================================================

class CustomAssertions(unittest.TestCase):

    def assert_status_of_phase(self, output, status, phase, test_name):
        """Asserts that 'output' contains a line showing the given
        status for the given phase for the given test_name"""
        expected = re.compile(r'^ *{} +'.format(re.escape(status)) +
                              self._test_name_and_phase_regex(test_name, phase),
                              flags=re.MULTILINE)

        six.assertRegex(self, output, expected)

    def assert_phase_absent(self, output, phase, test_name):
        """Asserts that 'output' does not contain a status line for the
        given phase and test_name"""
        expected = re.compile(r'^.* +' +
                              self._test_name_and_phase_regex(test_name, phase),
                              flags=re.MULTILINE)

        six_additions.assertNotRegex(self, output, expected)

    def assert_core_phases(self, output, test_name, fails):
        """Asserts that 'output' contains a line for each of the core test
        phases for the given test_name. All results should be PASS
        except those given by the fails list, which should be FAILS.
        """
        for phase in test_status.CORE_PHASES:
            if phase in fails:
                status = test_status.TEST_FAIL_STATUS
            else:
                status = test_status.TEST_PASS_STATUS
            self.assert_status_of_phase(output=output,
                                        status=status,
                                        phase=phase,
                                        test_name=test_name)

    @staticmethod
    def _test_name_and_phase_regex(test_name, phase):
        """Returns a regex matching the portion of a TestStatus line
        containing the test name and phase"""
        # The main purpose of extracting this into a shared method is:
        # assert_phase_absent could wrongly pass if the format of the
        # TestStatus output changed without that method's regex
        # changing. By making its regex shared as much as possible with
        # the regex in assert_status_of_phase, we decrease the chances
        # of these false passes.
        return r'{} +{} *$'.format(re.escape(test_name), re.escape(phase))

class TestCustomAssertions(CustomAssertions):

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

# ========================================================================
# Actual test class
# ========================================================================

class TestCsStatus(CustomAssertions):

    # ------------------------------------------------------------------------
    # Test helper functions
    # ------------------------------------------------------------------------

    # An arbitrary phase we can use when we want to work with a non-core phase
    _NON_CORE_PHASE = test_status.MEMLEAK_PHASE

    # Another arbitrary phase if we need two different non-core phases
    _NON_CORE_PHASE2 = test_status.BASELINE_PHASE

    def setUp(self):
        self._testroot = tempfile.mkdtemp()
        self._output = six.StringIO()

    def tearDown(self):
        self._output.close()
        shutil.rmtree(self._testroot, ignore_errors=True)

    def create_test_dir(self, test_dir):
        """Creates the given test directory under testroot.

        Returns the full path to the created test directory.
        """
        fullpath = os.path.join(self._testroot, test_dir)
        os.makedirs(fullpath)
        return fullpath

    @staticmethod
    def create_test_status_core_passes(test_dir_path, test_name):
        """Creates a TestStatus file in the given path, with PASS status
        for all core phases"""
        with test_status.TestStatus(test_dir=test_dir_path,
                                    test_name=test_name) as ts:
            for phase in test_status.CORE_PHASES:
                ts.set_status(phase, test_status.TEST_PASS_STATUS)

    def set_last_core_phase_to_fail(self, test_dir_path, test_name):
        """Sets the last core phase to FAIL

        Returns the name of this phase"""
        fail_phase = test_status.CORE_PHASES[-1]
        self.set_phase_to_status(test_dir_path=test_dir_path,
                                 test_name=test_name,
                                 phase=fail_phase,
                                 status=test_status.TEST_FAIL_STATUS)
        return fail_phase

    @staticmethod
    def set_phase_to_status(test_dir_path, test_name, phase, status):
        """Sets the given phase to the given status for this test"""
        with test_status.TestStatus(test_dir=test_dir_path,
                                    test_name=test_name) as ts:
            ts.set_status(phase, status)

    # ------------------------------------------------------------------------
    # Begin actual tests
    # ------------------------------------------------------------------------

    def test_single_test(self):
        """cs_status for a single test should include some minimal expected output"""
        test_name = 'my.test.name'
        test_dir = 'my.test.name.testid'
        test_dir_path = self.create_test_dir(test_dir)
        self.create_test_status_core_passes(test_dir_path, test_name)
        cs_status([os.path.join(test_dir_path, 'TestStatus')],
                  out=self._output)
        self.assert_core_phases(self._output.getvalue(), test_name, fails=[])

    def test_two_tests(self):
        """cs_status for two tests (one with a FAIL) should include some minimal expected output"""
        test_name1 = 'my.test.name1'
        test_name2 = 'my.test.name2'
        test_dir1 = test_name1 + '.testid'
        test_dir2 = test_name2 + '.testid'
        test_dir_path1 = self.create_test_dir(test_dir1)
        test_dir_path2 = self.create_test_dir(test_dir2)
        self.create_test_status_core_passes(test_dir_path1, test_name1)
        self.create_test_status_core_passes(test_dir_path2, test_name2)
        test2_fail_phase = self.set_last_core_phase_to_fail(test_dir_path2, test_name2)
        cs_status([os.path.join(test_dir_path1, 'TestStatus'),
                   os.path.join(test_dir_path2, 'TestStatus')],
                  out=self._output)
        self.assert_core_phases(self._output.getvalue(), test_name1, fails=[])
        self.assert_core_phases(self._output.getvalue(), test_name2, fails=[test2_fail_phase])

    def test_fails_only(self):
        """With fails_only flag, only fails and pends should appear in the output"""
        test_name = 'my.test.name'
        test_dir = 'my.test.name.testid'
        test_dir_path = self.create_test_dir(test_dir)
        self.create_test_status_core_passes(test_dir_path, test_name)
        fail_phase = self.set_last_core_phase_to_fail(test_dir_path, test_name)
        pend_phase = self._NON_CORE_PHASE
        self.set_phase_to_status(test_dir_path, test_name,
                                 phase=pend_phase,
                                 status=test_status.TEST_PEND_STATUS)
        cs_status([os.path.join(test_dir_path, 'TestStatus')],
                  fails_only=True,
                  out=self._output)
        self.assert_status_of_phase(output=self._output.getvalue(),
                                    status=test_status.TEST_FAIL_STATUS,
                                    phase=fail_phase,
                                    test_name=test_name)
        self.assert_status_of_phase(output=self._output.getvalue(),
                                    status=test_status.TEST_PEND_STATUS,
                                    phase=pend_phase,
                                    test_name=test_name)
        for phase in test_status.CORE_PHASES:
            if phase != fail_phase:
                self.assert_phase_absent(output=self._output.getvalue(),
                                         phase=phase,
                                         test_name=test_name)
        six_additions.assertNotRegex(self, self._output.getvalue(), r'Overall:')

    def test_count_fails(self):
        """Test the count of fails with three tests

        For first phase of interest: First test FAILs, second PASSes,
        third FAILs; count should be 2, and this phase should not appear
        individually for each test.

        For second phase of interest: First test PASSes, second PASSes,
        third FAILs; count should be 1, and this phase should not appear
        individually for each test.
        """
        # Note that this test does NOT cover:
        # - combining count_fails_phase_list with fails_only: currently,
        #   this wouldn't cover any additional code/logic
        # - ensuring that PENDs are also counted: currently, this
        #   wouldn't cover any additional code/logic
        phase_of_interest1 = self._NON_CORE_PHASE
        phase_of_interest2 = self._NON_CORE_PHASE2
        statuses1 = [test_status.TEST_FAIL_STATUS,
                     test_status.TEST_PASS_STATUS,
                     test_status.TEST_FAIL_STATUS]
        statuses2 = [test_status.TEST_PASS_STATUS,
                     test_status.TEST_PASS_STATUS,
                     test_status.TEST_FAIL_STATUS]
        test_paths = []
        test_names = []
        for testnum in range(3):
            test_name = 'my.test.name' + str(testnum)
            test_names.append(test_name)
            test_dir = test_name + '.testid'
            test_dir_path = self.create_test_dir(test_dir)
            self.create_test_status_core_passes(test_dir_path, test_name)
            self.set_phase_to_status(test_dir_path, test_name,
                                     phase=phase_of_interest1,
                                     status=statuses1[testnum])
            self.set_phase_to_status(test_dir_path, test_name,
                                     phase=phase_of_interest2,
                                     status=statuses2[testnum])
            test_paths.append(os.path.join(test_dir_path, 'TestStatus'))

        cs_status(test_paths,
                  count_fails_phase_list=[phase_of_interest1, phase_of_interest2],
                  out=self._output)

        for testnum in range(3):
            self.assert_phase_absent(output=self._output.getvalue(),
                                     phase=phase_of_interest1,
                                     test_name=test_names[testnum])
            self.assert_phase_absent(output=self._output.getvalue(),
                                     phase=phase_of_interest2,
                                     test_name=test_names[testnum])
        count_regex1 = r'{} +non-passes: +2'.format(re.escape(phase_of_interest1))
        six.assertRegex(self, self._output.getvalue(), count_regex1)
        count_regex2 = r'{} +non-passes: +1'.format(re.escape(phase_of_interest2))
        six.assertRegex(self, self._output.getvalue(), count_regex2)

if __name__ == '__main__':
    unittest.main()
