"""
This module contains a class that extends unittest.TestCase, adding custom assertions that
can be used when testing TestStatus.
"""

from CIME.XML.standard_module_setup import *

import unittest
import re
import six
import six_additions
from CIME import test_status

class CustomAssertionsTestStatus(unittest.TestCase):

    def assert_status_of_phase(self, output, status, phase, test_name, xfail=None):
        """Asserts that 'output' contains a line showing the given
        status for the given phase for the given test_name.

        'xfail' should have one of the following values:
        - None (the default): assertion passes regardless of whether there is an
          EXPECTED/UNEXPECTED string
        - 'no': The line should end with the phase, with no additional text after that
        - 'expected': After the phase, the line should contain '(EXPECTED FAILURE)'
        - 'unexpected': After the phase, the line should contain '(UNEXPECTED'
        """
        expected = (r'^ *{} +'.format(re.escape(status)) +
                    self._test_name_and_phase_regex(test_name, phase))

        if xfail == 'no':
            # There should be no other text after the testname and phase regex
            expected += r' *$'
        elif xfail == 'expected':
            expected += r' *{}'.format(re.escape(test_status.TEST_EXPECTED_FAILURE_COMMENT))
        elif xfail == 'unexpected':
            expected += r' *{}'.format(re.escape(test_status.TEST_UNEXPECTED_FAILURE_COMMENT_START))
        else:
            expect(xfail is None, "Unhandled value of xfail argument")

        expected_re = re.compile(expected, flags=re.MULTILINE)
        six.assertRegex(self, output, expected_re)

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

    def assert_num_expected_unexpected_fails(self, output, num_expected, num_unexpected):
        """Asserts that the number of occurrences of expected and unexpected fails in
        'output' matches the given numbers"""
        self.assertEqual(output.count(test_status.TEST_EXPECTED_FAILURE_COMMENT), num_expected)
        self.assertEqual(output.count(test_status.TEST_UNEXPECTED_FAILURE_COMMENT_START), num_unexpected)

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
        return r'{} +{}'.format(re.escape(test_name), re.escape(phase))
