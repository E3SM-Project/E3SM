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


from rtest_betr import Tolerances


class Tolerances_suite(unittest.TestCase):
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

    def test_tolerances_default(self):
        """Test basic initialization of a tolerance object and return of a
        default value.

        """
        tolerances = Tolerances()
        expected = Tolerances._DEFAULT_EPSILON
        received = tolerances.get(Tolerances.CONC, 'value')

        self.assertEqual(expected, received)

    # ------------------------------------------------------

    def test_tolerances_check_valid_category_valid(self):
        """Check that a user specified category is valid.
        """
        tolerances = Tolerances()
        category = Tolerances.CONC

        expected = tolerances.check_valid_category(category)
        self.assertTrue(expected)

    def test_tolerances_check_valid_category_invalid(self):
        """Check that an error is raised for an invalid user specified
        category

        """
        tolerances = Tolerances()
        category = 'junk'

        self.assertRaises(RuntimeError,
                          tolerances.check_valid_category, category)

    # ------------------------------------------------------

    def test_tolerances_update_name_valid(self):
        """Test that the update method accepts valid data.

        """
        tolerances = Tolerances()
        category = Tolerances.GENERAL
        expected_value = '1.0e-2'
        expected_type = 'absolute'
        data = ' '.join([expected_value, expected_type])
        tolerances.update_from_name(category, data)

        # NOTE(bja, 201603) bad test, depends on object internals.
        received_value = tolerances._tolerances[category]['value']
        received_type = tolerances._tolerances[category]['type']
        self.assertEqual(float(expected_value), received_value)
        self.assertEqual(expected_type, received_type)

    def test_tolerances_update_name_invalid_type(self):
        """Test that the update method raises an error for an unknown type.

        """
        tolerances = Tolerances()
        category = 'foo'
        expected_value = '1.0e-3'
        expected_type = 'absolute'
        data = ' '.join([expected_value, expected_type])
        self.assertRaises(RuntimeError,
                          tolerances.update_from_name, category, data)

    def test_tolerances_update_name_invalid_float(self):
        """Test that the update method raises an error for a value that can't
be converted to a float.

        """
        tolerances = Tolerances()
        category = Tolerances.GENERAL
        expected_value = 'one.two'
        expected_type = 'absolute'
        data = ' '.join([expected_value, expected_type])
        self.assertRaises(RuntimeError,
                          tolerances.update_from_name, category, data)

    def test_tolerances_update_name_invalid_minimum(self):
        """Test that the update from name method raises an error for a value
less than zero.

        """
        tolerances = Tolerances()
        category = Tolerances.GENERAL
        expected_value = '-1.0e-2'
        expected_type = 'absolute'
        data = ' '.join([expected_value, expected_type])
        self.assertRaises(RuntimeError,
                          tolerances.update_from_name, category, data)

    def test_tolerances_update_name_invalid_max(self):
        """Test that the update from name method raises an error for a value
that exceeds the maximum.

        """
        tolerances = Tolerances()
        category = Tolerances.GENERAL
        expected_value = '1.0e20000'
        expected_type = 'absolute'
        data = ' '.join([expected_value, expected_type])
        self.assertRaises(RuntimeError,
                          tolerances.update_from_name, category, data)

    def test_tolerances_update_name_invalid_max_percent(self):
        """Test that the update from name raises an error for a percent > 100

        """
        tolerances = Tolerances()
        category = Tolerances.GENERAL
        expected_value = '100.1'
        expected_type = 'percent'
        data = ' '.join([expected_value, expected_type])
        self.assertRaises(RuntimeError,
                          tolerances.update_from_name, category, data)

    # ------------------------------------------------------

    def test_tolerances_get_valid(self):
        """Test that the get method return the correct value for a requested key

        """
        tolerances = Tolerances()
        category = Tolerances.CONC
        expected_value = '1.2345'
        expected_type = 'percent'
        data = ' '.join([expected_value, expected_type])
        tolerances.update_from_name(category, data)
        received_value = tolerances.get(category, 'value')
        received_type = tolerances.get(category, 'type')
        received_max = tolerances.get(category, 'max')
        self.assertEqual(float(expected_value), received_value)
        self.assertEqual(expected_type, received_type)
        self.assertEqual(100.0, received_max)

    def test_tolerances_get_invalid_category(self):
        """Test that the get method raises an error for a invalid category.

        """
        tolerances = Tolerances()
        category = 'junk'
        key = 'min'
        # have to manually set the internals because error checking in
        # the utility functions catches this error too early.
        tolerances._tolerances[category] = {}
        self.assertRaises(RuntimeError,
                          tolerances.get, category, key)

    def test_tolerances_get_invalid_key(self):
        """Test that the get method raises an error for an invalid tolerance
data key.  type.

        """
        tolerances = Tolerances()
        category = Tolerances.GENERAL
        key = 'junk'
        self.assertRaises(RuntimeError,
                          tolerances.get, category, key)


if __name__ == '__main__':
    # unittest.main(buffer=True)
    unittest.main()  # pragma: no coverage
