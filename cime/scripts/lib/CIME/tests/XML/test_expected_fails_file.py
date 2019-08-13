#!/usr/bin/env python

import unittest
import os
import shutil
import tempfile
import six
from CIME.XML.expected_fails_file import ExpectedFailsFile
from CIME.utils import CIMEError
from CIME.expected_fails import ExpectedFails

class TestExpectedFailsFile(unittest.TestCase):

    def setUp(self):
        self._workdir = tempfile.mkdtemp()
        self._xml_filepath = os.path.join(self._workdir, "expected_fails.xml")

    def tearDown(self):
        shutil.rmtree(self._workdir)

    def test_basic(self):
        """Basic test of the parsing of an expected fails file"""
        contents = """<?xml version= "1.0"?>
<expectedFails version="1.1">
  <test name="my.test.1">
    <phase name="RUN">
      <status>FAIL</status>
      <issue>#404</issue>
    </phase>
    <phase name="COMPARE_base_rest">
      <status>PEND</status>
      <issue>#404</issue>
      <comment>Because of the RUN failure, this phase is listed as PEND</comment>
    </phase>
  </test>
  <test name="my.test.2">
    <phase name="GENERATE">
      <status>FAIL</status>
      <issue>ESMCI/cime#2917</issue>
    </phase>
    <phase name="BASELINE">
      <status>FAIL</status>
      <issue>ESMCI/cime#2917</issue>
    </phase>
  </test>
</expectedFails>
"""
        with open(self._xml_filepath, 'w') as xml_file:
            xml_file.write(contents)
        expected_fails_file = ExpectedFailsFile(self._xml_filepath)
        xfails = expected_fails_file.get_expected_fails()

        expected_test1 = ExpectedFails()
        expected_test1.add_failure('RUN', 'FAIL')
        expected_test1.add_failure('COMPARE_base_rest', 'PEND')
        expected_test2 = ExpectedFails()
        expected_test2.add_failure('GENERATE', 'FAIL')
        expected_test2.add_failure('BASELINE', 'FAIL')
        expected = {'my.test.1': expected_test1,
                    'my.test.2': expected_test2}

        self.assertEqual(xfails, expected)

    def test_same_test_appears_twice(self):
        """If the same test appears twice, its information should be appended.

        This is not the typical, expected layout of the file, but it should be handled
        correctly in case the file is written this way.
        """
        contents = """<?xml version= "1.0"?>
<expectedFails version="1.1">
  <test name="my.test.1">
    <phase name="RUN">
      <status>FAIL</status>
      <issue>#404</issue>
    </phase>
  </test>
  <test name="my.test.1">
    <phase name="COMPARE_base_rest">
      <status>PEND</status>
      <issue>#404</issue>
      <comment>Because of the RUN failure, this phase is listed as PEND</comment>
    </phase>
  </test>
</expectedFails>
"""
        with open(self._xml_filepath, 'w') as xml_file:
            xml_file.write(contents)
        expected_fails_file = ExpectedFailsFile(self._xml_filepath)
        xfails = expected_fails_file.get_expected_fails()

        expected_test1 = ExpectedFails()
        expected_test1.add_failure('RUN', 'FAIL')
        expected_test1.add_failure('COMPARE_base_rest', 'PEND')
        expected = {'my.test.1': expected_test1}

        self.assertEqual(xfails, expected)

    def test_invalid_file(self):
        """Given an invalid file, an exception should be raised in schema validation"""

        # This file is missing a <status> element in the <phase> block.
        #
        # It's important to have the expectedFails version number be greater than 1,
        # because schema validation isn't done in cime for files with a version of 1.
        contents = """<?xml version= "1.0"?>
<expectedFails version="1.1">
  <test name="my.test.1">
    <phase name="RUN">
    </phase>
  </test>
</expectedFails>
"""
        with open(self._xml_filepath, 'w') as xml_file:
            xml_file.write(contents)

        with six.assertRaisesRegex(self, CIMEError, "Schemas validity error"):
            _ = ExpectedFailsFile(self._xml_filepath)

if __name__ == '__main__':
    unittest.main()
