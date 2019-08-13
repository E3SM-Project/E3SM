#!/usr/bin/env python

import unittest
from CIME.utils import indent_string

class TestIndentStr(unittest.TestCase):
    """Test the indent_string function.

    """

    def test_indent_string_singleline(self):
        """Test the indent_string function with a single-line string

        """
        mystr = 'foo'
        result = indent_string(mystr, 4)
        expected = '    foo'
        self.assertEqual(expected, result)

    def test_indent_string_multiline(self):
        """Test the indent_string function with a multi-line string

        """
        mystr = """hello
hi
goodbye
"""
        result = indent_string(mystr, 2)
        expected = """  hello
  hi
  goodbye
"""
        self.assertEqual(expected, result)

if __name__ == '__main__':
    unittest.main()

