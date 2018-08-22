#!/usr/bin/env python

import unittest
import shutil
import os
import tempfile
import StringIO
import re
import six
from CIME.cs_status import cs_status
from CIME import test_status

# ========================================================================
# Custom assertions
# ========================================================================

class CustomAssertions(unittest.TestCase):

    def assert_core_phases(self, output, test_name, fails):
        """Asserts that 'output' contains a line for each of the core test
        phases for the given test_name. All results should be PASS
        except those given by the fails list, which should be FAILS.
        """
        for phase in test_status.CORE_PHASES:
            if phase in fails:
                status = 'FAIL'
            else:
                status = 'PASS'
            expected = r'{} +{} +{}'.format(status, re.escape(test_name), re.escape(phase))
            six.assertRegex(self, output, expected)

class TestCustomAssertions(CustomAssertions):

    def output_line(self, status, test_name, phase):
        return status + ' ' + test_name + ' ' + phase + '\n'

    def test_assert_core_phases_passes(self):
        """assert_core_phases passes when it should"""
        output = ''
        fails = [test_status.CORE_PHASES[1]]
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            if phase in fails:
                status = 'FAIL'
            else:
                status = 'PASS'
            output = output + self.output_line(status, test_name, phase)

        self.assert_core_phases(output, test_name, fails)

    def test_assert_core_phases_missing_phase_fails(self):
        """assert_core_phases fails if there is a missing phase"""
        output = ''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            if phase != test_status.CORE_PHASES[1]:
                output = output + self.output_line('PASS', test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, test_name, fails=[])

    def test_assert_core_phases_wrong_status_fails(self):
        """assert_core_phases fails if a phase has the wrong status"""
        output =''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            output = output + self.output_line('PASS', test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, test_name, fails=[test_status.CORE_PHASES[1]])

    def test_assert_core_phases_wrong_name_fails(self):
        """assert_core_phases fails if the test name is wrong"""
        output =''
        test_name = 'my.test.name'
        for phase in test_status.CORE_PHASES:
            output = output + self.output_line('PASS', test_name, phase)

        with self.assertRaises(AssertionError):
            self.assert_core_phases(output, 'my.test', fails=[])

# ========================================================================
# Actual test class
# ========================================================================

class TestCsStatus(CustomAssertions):

    # ------------------------------------------------------------------------
    # Test helper functions
    # ------------------------------------------------------------------------

    def setUp(self):
        self._testroot = tempfile.mkdtemp()
        self._output = StringIO.StringIO()

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

    def create_test_status_core_passes(self, test_dir_path, test_name):
        """Creates a TestStatus file in the given path, with PASS status
        for all core phases"""
        with test_status.TestStatus(test_dir=test_dir_path,
                                    test_name=test_name) as ts:
            for phase in test_status.CORE_PHASES:
                ts.set_status(phase, test_status.TEST_PASS_STATUS)

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
        test2_fail_phase = test_status.CORE_PHASES[-1]
        self.create_test_status_core_passes(test_dir_path1, test_name1)
        self.create_test_status_core_passes(test_dir_path2, test_name2)
        with test_status.TestStatus(test_dir=test_dir_path2,
                                    test_name=test_name2) as ts:
            ts.set_status(test2_fail_phase, test_status.TEST_FAIL_STATUS)
        cs_status([os.path.join(test_dir_path1, 'TestStatus'),
                   os.path.join(test_dir_path2, 'TestStatus')],
                  out=self._output)
        self.assert_core_phases(self._output.getvalue(), test_name1, fails=[])
        self.assert_core_phases(self._output.getvalue(), test_name2, fails=[test2_fail_phase])

if __name__ == '__main__':
    unittest.main()
