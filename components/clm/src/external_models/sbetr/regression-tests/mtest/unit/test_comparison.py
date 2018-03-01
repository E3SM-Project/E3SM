#!/usr/bin/env python
"""
Unit test suite for the betr regression test manager

"""
from __future__ import print_function

import logging
import os
import os.path
import sys
import unittest

if sys.version_info[0] == 2:  # pragma: no cover
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser


from rtest_betr import Comparison


class Comparison_suite(unittest.TestCase):
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

    def test_comparison_get_section_category_valid(self):
        """Test that we can extract the category correctly from a test
        description section.

        """
        conf = {}
        comparison = Comparison('unittest', conf)

        section = {'category': comparison._tolerances.CONC,
                   'foo': 'bar baz'}
        expected = comparison._tolerances.CONC
        received = comparison._get_section_category('some test', section)
        self.assertEqual(expected, received)

    def test_comparison_get_section_category_missing_category(self):
        """Test that we can extract the category correctly from a test
        description section.

        """
        conf = {}
        comparison = Comparison('unittest', conf)

        section = {'pet': 'dog cat',
                   'foo': 'bar baz'}

        self.assertRaises(RuntimeError,
                          comparison._get_section_category,
                          'some test', section)

    def test_comparison_get_section_category_invalid_category(self):
        """Test that we can extract the category correctly from a test
        description section.

        """
        conf = {}
        comparison = Comparison('unittest', conf)

        section = {'category': 'dog cat',
                   'foo': 'bar baz'}

        self.assertRaises(RuntimeError,
                          comparison._get_section_category,
                          'some test', section)

    # ------------------------------------------------------
    def test_comparison_float_absolute_pass(self):
        """Test that comparison of float with absolute tolerance passes when
        it is less that tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '1.0e-16 absolute')
        section = 'Ca'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.1e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertTrue(received)

    def test_comparison_float_absolute_fail(self):
        """Test that comparison of float with absolute tolerance fails when it
        is greater than tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '1.0e-18 absolute')
        section = 'Ca'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.1e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertFalse(received)

    # ------------------------------------------------------
    def test_comparison_float_relative_pass(self):
        """Test that comparison of float with relative tolerance passes when
        it is less that tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '1.0e-4 relative')
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.00001e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertTrue(received)

    def test_comparison_float_relative_fail(self):
        """Test that comparison of float with relative tolerance fails when it
        is greater than tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '1.0e-5 relative')
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.00001e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertFalse(received)

    # ------------------------------------------------------
    def test_comparison_float_percent_pass(self):
        """Test that comparison of float with percent tolerance passes when
        it is less that tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '5.0 percent')
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.04e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertTrue(received)

    def test_comparison_float_percent_fail(self):
        """Test that comparison of float with percent tolerance fails when it
        is greater than tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '5.0 percent')
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.06e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertFalse(received)

    # ------------------------------------------------------
    def test_comparison_float_invalid_method(self):
        """Test that comparison of float with junk tolerance type fails with
        an exception. Kind of a contrived example because it should be
        very hard to get to this code given all the prior error
        checking.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '5.0 percent')
        comparison._tolerances._tolerances[category]['type'] = 'junk'
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.06e-16'

        self.assertRaises(RuntimeError,
                          comparison._compare_float_values_with_tolerance,
                          category, section, key, a_data, b_data)

    def test_comparison_integer_unimplemented(self):
        """Test that comparison of integer/discrete fails with with runtime
        error because it hasn't been implemented yet.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = 'discrete'
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.06e-16'

        self.assertRaises(RuntimeError,
                          comparison._compare_integer_values_with_tolerance,
                          category, section, key, a_data, b_data)

    # ------------------------------------------------------
    def test_compare_values_pass(self):
        """Test that call to compare_values_with_tolerance with valid data
        passes.

        FIXME(bja, 201603) Only way to check this is via the _status
        flag being unset... very bad.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '1.0e-4 relative')
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.00001e-16'

        comparison._compare_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertIsNone(comparison._status)

    def test_compare_values_fail(self):
        """Test that call to compare_values_with_tolerance fails with
        unimplemented integer data

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = 'discrete'
        section = 'foo'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.06e-16'

        self.assertRaises(RuntimeError,
                          comparison._compare_values_with_tolerance,
                          category, section, key, a_data, b_data)

    # ------------------------------------------------------
    def test_compare_options_different_categories(self):
        """Test that a regression sections having different data categories
        causes a failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        section = 'foo_test'
        a_data = {
            'category': comparison._tolerances.GENERAL,
            'min': '1.02345e-1',
        }
        a_name = 'a_baseline'
        b_data = {
            'category': comparison._tolerances.CONC,
            'min': '1.02345e-1',
        }
        b_name = 'b_regression'
        comparison._compare_options(section, a_data, a_name, b_data, b_name)
        self.assertEqual(comparison._status, 'fail')

    def test_compare_options_missing_key_first_fails(self):
        """Test that a key missing from the first regression section causes a
        failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        section = 'foo_test'
        a_data = {
            'category': comparison._tolerances.CONC,
            'min': '1.02345e-1',
            'max': '9.8765e1',
        }
        a_name = 'a_baseline'
        b_data = {
            'category': comparison._tolerances.CONC,
            'min': '1.02345e-1',
        }
        b_name = 'b_regression'
        comparison._compare_options(section, a_data, a_name, b_data, b_name)
        self.assertEqual(comparison._status, 'fail')

    def test_compare_options_missing_key_second_fails(self):
        """Test that a key missing from the second regression section causes a
        failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        section = 'foo_test'
        a_data = {
            'category': comparison._tolerances.CONC,
            'min': '1.02345e-1',
        }
        a_name = 'a_baseline'
        b_data = {
            'category': comparison._tolerances.CONC,
            'min': '1.02345e-1',
            'max': '9.8765e1',
        }
        b_name = 'b_regression'
        comparison._compare_options(section, a_data, a_name, b_data, b_name)
        self.assertEqual(comparison._status, 'fail')

    # ------------------------------------------------------
    def test_compare_update_valid(self):
        """Test that the comparison object passes tolerance update info
        through to the tolerances object.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        value = '123.4'
        category = comparison._tolerances.CONC
        new_conf = {'default_tolerances':
                    {
                        category: '123.4 relative',
                    },
        }
        comparison.update(new_conf)
        self.assertEqual(comparison._tolerances._tolerances[category]['value'],
                         float(value))

    # ------------------------------------------------------
    def test_set_denominator_first_zero(self):
        """Test that setting the denominator works when the first value is
        zero.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = 0.0
        b_data = 1.00001e-16

        denominator = comparison._set_denominator(a_data, b_data)
        self.assertEqual(b_data, denominator)

    def test_set_denominator_second_zero(self):
        """Test that setting the denominator works when the second value is
        zero.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = 1.00001e-16
        b_data = 0.0

        denominator = comparison._set_denominator(a_data, b_data)
        self.assertEqual(a_data, denominator)

    def test_set_denominator_both_zero(self):
        """Test that setting the denominator works when the both values are
        zero.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = 0.0
        b_data = 0.0

        denominator = comparison._set_denominator(a_data, b_data)
        self.assertEqual(1.0, denominator)

    def test_set_denominator_neither_zero(self):
        """Test that setting the denominator returns the first value when
        the both values are non-zero.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = 1.0e-1
        b_data = 2.0e-12

        denominator = comparison._set_denominator(a_data, b_data)
        self.assertEqual(a_data, denominator)

    # ------------------------------------------------------
    def test_compare_sections_missing_baseline_fails(self):
        """Test that a key missing from baseline regression section causes a
        failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        section = 'foo_test'
        a_data = {
            section: {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
            },
        }
        a_name = 'a_baseline'
        b_data = {
            section: {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {},
        }
        b_name = 'b_regression'
        comparison._compare_sections(a_data, a_name, b_data, b_name)
        self.assertEqual(comparison._status, 'fail')

    def test_compare_sections_missing_current_fails(self):
        """Test that a key missing from current regression section causes a
        failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        section = 'foo_test'
        a_data = {
            section: {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
            },
            'other_section': {},
        }
        a_name = 'a_baseline'
        b_data = {
            section: {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
        }
        b_name = 'b_regression'
        comparison._compare_sections(a_data, a_name, b_data, b_name)
        self.assertEqual(comparison._status, 'fail')

    def test_compare_sections_both_pass(self):
        """Test that if all sections are in both regression and baseline
        doesn't cause a failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        section = 'foo_test'
        a_data = {
            section: {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }
        a_name = 'a_baseline'
        b_data = {
            section: {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }
        b_name = 'b_regression'
        comparison._compare_sections(a_data, a_name, b_data, b_name)
        self.assertIsNone(comparison._status)

    # ------------------------------------------------------

    def test_regression_to_baseline_pass(self):
        """Test that if regression and baseline are the same, we get a pass

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = {
            'section': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }
        b_data = {
            'section': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }

        status = comparison.regression_to_baseline(a_data, b_data)
        self.assertEqual(status, 'pass')

    def test_regression_to_baseline_fail_tolerance(self):
        """Test that if regression and baseline have values that exceed
        tolerances, we get a failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = {
            'section': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }
        b_data = {
            'section': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '3.45e-2',
            },
        }

        status = comparison.regression_to_baseline(a_data, b_data)
        self.assertEqual(status, 'fail')

    def test_regression_to_baseline_fail_key(self):
        """Test that if regression and baseline don't have the same keys, we
        get a failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = {
            'section': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'max': '2.345e-2',
            },
        }
        b_data = {
            'section': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }

        status = comparison.regression_to_baseline(a_data, b_data)
        self.assertEqual(status, 'fail')

    def test_regression_to_baseline_fail_section(self):
        """Test that if regression and baseline don't have the same sections,
        we get a failure.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        a_data = {
            'foo': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }
        b_data = {
            'section': {
                'category': comparison._tolerances.CONC,
                'min': '1.02345e-1',
                'max': '9.8765e1',
            },
            'other_section': {
                'category': comparison._tolerances.CONC,
                'min': '2.345e-2',
            },
        }

        status = comparison.regression_to_baseline(a_data, b_data)
        self.assertEqual(status, 'fail')

    # ------------------------------------------------------

if __name__ == '__main__':
    # unittest.main(buffer=True)
    unittest.main()  # pragma: no cover
