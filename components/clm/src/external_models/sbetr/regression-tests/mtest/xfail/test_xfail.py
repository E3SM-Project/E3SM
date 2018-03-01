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

from rtest_betr import main, commandline_options


class BeTR_suite(unittest.TestCase):
    """
    """
    _LOG_FILENAME = 'dummy.testlog'

    def setUp(self):
        """
        """
        test_dir = 'mtest/xfail'
        ext_args = [
            '--executable', '{0}/dummy.exe'.format(test_dir),
            '--config', '{0}/xfail.cfg'.format(test_dir),
            '--check-only',
        ]
        self._options = commandline_options(ext_args)

    def tearDown(self):
        """
        """
        pass

    # ------------------------------------------------------

    def test_example_default(self):
        """A working example test
        default value.

        """
        status = main(self._options)
        self.assertNotEqual(status, 0)

if __name__ == '__main__':
    # unittest.main(buffer=True)
    unittest.main()  # pragma: no coverage
