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

from rtest_betr import main


class BeTR_suite(unittest.TestCase):
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

    def test_example_default(self):
        """A working example test
        default value.

        """
        expected = True
        received = True

        self.assertEqual(expected, received)


if __name__ == '__main__':
    # unittest.main(buffer=True)
    unittest.main()  # pragma: no coverage
