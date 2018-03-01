#!/usr/bin/env python
"""
Unit test suite for the betr regression test manager

"""
from __future__ import print_function

import logging
import os
import sys
import unittest

if sys.version_info[0] == 2:  # pragma: no coverage
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser


from rtest_betr import StatusSummary


class StatusSummary_suite(unittest.TestCase):
    """
    """
    _LOG_FILENAME = 'dummy.testlog'

    def setUp(self):
        """
        """
        logging.basicConfig(filename=self._LOG_FILENAME,
                            filemode='w',
                            level=logging.INFO,
                            format='%(message)s')
        logging.info('mtest {0} unit test log.'.format(__name__))

    def tearDown(self):
        """
        """
        logging.shutdown()
        if os.path.isfile(self._LOG_FILENAME):
            os.remove(self._LOG_FILENAME)  # pragma: no coverage

    # ------------------------------------------------------

    def test_status_default(self):
        """Test basic initialization of a status summar object and return of a
        default value.

        """
        status = StatusSummary()
        expected = 0

        received = status.skips()
        self.assertEqual(expected, received)

        received = status.failures()
        self.assertEqual(expected, received)

        received = status.passes()
        self.assertEqual(expected, received)

        received = status.total()
        self.assertEqual(expected, received)

    def test_status_add_skip_default(self):
        """Test default increment of skips

        """
        status = StatusSummary()
        expected = 2

        status.add_skip()
        status.add_skip()
        received = status.skips()
        self.assertEqual(expected, received)

    def test_status_add_skip_specified(self):
        """Test specified increment of skips

        """
        status = StatusSummary()

        expected = 2
        status.add_skip(2)
        received = status.skips()
        self.assertEqual(expected, received)

        expected = 5
        status.add_skip(3)
        received = status.skips()
        self.assertEqual(expected, received)

    def test_status_add_fail_default(self):
        """Test default increment of failures

        """
        status = StatusSummary()
        expected = 2

        status.add_failure()
        status.add_failure()
        received = status.failures()
        self.assertEqual(expected, received)

    def test_status_add_failures_specified(self):
        """Test default increment of failures

        """
        status = StatusSummary()

        expected = 2
        status.add_failure(2)
        received = status.failures()
        self.assertEqual(expected, received)

        expected = 5
        status.add_failure(3)
        received = status.failures()
        self.assertEqual(expected, received)

    def test_status_add_pass_default(self):
        """Test default increment of passes

        """
        status = StatusSummary()
        expected = 3

        status.add_pass()
        status.add_pass()
        status.add_pass()
        received = status.passes()
        self.assertEqual(expected, received)

    def test_status_add_pass_specified(self):
        """Test default increment of passs

        """
        status = StatusSummary()

        expected = 2
        status.add_pass(2)
        received = status.passes()
        self.assertEqual(expected, received)

        expected = 6
        status.add_pass(4)
        received = status.passes()
        self.assertEqual(expected, received)

    def test_status_add_total_default(self):
        """Test default increment of total when passes, failures and skips are
added.

        """
        status = StatusSummary()
        expected = 3

        status.add_pass()
        status.add_failure()
        status.add_skip()
        received = status.total()
        self.assertEqual(expected, received)

    def test_status_add_total_specified(self):
        """Test increment of total when passs, failures and skips are added.

        """
        status = StatusSummary()

        expected = 2
        status.add_failure(2)
        received = status.total()
        self.assertEqual(expected, received)

        expected = 6
        status.add_pass(4)
        received = status.total()
        self.assertEqual(expected, received)

        expected = 9
        status.add_skip(3)
        received = status.total()
        self.assertEqual(expected, received)

    def test_status_addition_operator(self):
        """Test addition operator for adding status objects

        """
        status_1 = StatusSummary()
        expected_total_1 = 3

        status_1.add_pass()
        status_1.add_failure()
        status_1.add_skip()
        received = status_1.total()
        self.assertEqual(expected_total_1, received)

        status_2 = StatusSummary()

        status_2.add_failure(2)
        status_2.add_pass(4)
        status_2.add_skip(3)

        expected_total_2 = 9
        received = status_2.total()
        self.assertEqual(expected_total_2, received)

        status_3 = status_1 + status_2
        expected_total_3 = expected_total_1 + expected_total_2
        received = status_3.total()
        self.assertEqual(expected_total_3, received)

        expected_skip_3 = status_1.skips() + status_2.skips()
        received = status_3.skips()
        self.assertEqual(expected_skip_3, received)

    def test_status_addition_increment(self):
        """Test addition increment, +=, operator for adding status objects

        """
        status_1 = StatusSummary()
        expected_total_1 = 3

        status_1.add_pass()
        status_1.add_failure()
        status_1.add_skip()
        received = status_1.total()
        self.assertEqual(expected_total_1, received)

        status_2 = StatusSummary()

        status_2.add_failure(2)
        status_2.add_pass(4)
        status_2.add_skip(3)

        expected_total_2 = 9
        received = status_2.total()
        self.assertEqual(expected_total_2, received)

        expected_total_3 = expected_total_1 + expected_total_2
        expected_skip_3 = status_1.skips() + status_2.skips()

        status_1 += status_2

        received_skip_3 = status_1.skips()
        self.assertEqual(expected_skip_3, received_skip_3)

        received_total_3 = status_1.total()
        self.assertEqual(expected_total_3, received_total_3)

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    # unittest.main(buffer=True)
    unittest.main()  # pragma: no coverage
