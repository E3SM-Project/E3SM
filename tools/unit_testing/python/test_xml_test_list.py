#!/usr/bin/env python
"""Unit tests for the xml_test_list module.

Public classes:
TestTestSuiteSpec - TestSuiteSpec class unit tests.
TestSuitesFromXML - suites_from_xml unit tests.
"""

from os.path import abspath
import unittest
from xml.etree.ElementTree import XML, ElementTree

from xml_test_list import TestSuiteSpec, suites_from_xml

__all__ = ("TestSuitesFromXML", "TestSuitesFromXML")

class TestTestSuiteSpec(unittest.TestCase):

    """Tests for the TestSuiteSpec class."""

    def test_absolute_path(self):
        """TestSuiteSpec works as intended on absolute paths."""
        spec = TestSuiteSpec("name", [None], ["/path"])
        self.assertEqual("name", spec.name)
        self.assertEqual(["/path"], spec.directories)

    def test_relative_path(self):
        """TestSuiteSpec works as intended on relative paths."""
        spec = TestSuiteSpec("name", [None, None], ["path", "./path"])
        self.assertEqual([abspath("path"), abspath("./path")],
                         spec.directories)

    def test_no_path(self):
        """TestSuiteSpec works with no paths."""
        spec = TestSuiteSpec("name", [], [])
        self.assertEqual([], spec.directories)

    def test_label(self):
        """TestSuiteSpec takes and stores labels for input paths."""
        spec = TestSuiteSpec("name", ["foo", "bar"], ["/foo", "/bar"])
        self.assertEqual(["foo", "bar"], spec.labels)
        self.assertEqual(["/foo", "/bar"], spec.directories)

    def test_no_label(self):
        """TestSuiteSpec makes up a label for unlabeled input paths."""
        spec = TestSuiteSpec("name", [None, None], ["/foo", "/bar"])
        for label in spec.labels:
            self.assertEqual(label, TestSuiteSpec.UNLABELED_STRING)

    def test_iterate(self):
        """TestSuiteSpec provides an iterator over directories."""
        spec = TestSuiteSpec("name", ["foo", "bar"], ["/foo", "/bar"])
        self.assertEqual([("foo", "/foo"), ("bar", "/bar")],
                         list(d for d in spec))

class TestSuitesFromXML(unittest.TestCase):

    """Tests for the suites_from_xml function."""

    def check_spec_list(self, xml_str, names, directories,
                        known_paths={}, labels=None):
        """Check that a spec list matches input names and directories.

        This is used by the following tests to do the dirty work of making
        the list and making assertions about the names/paths.
        """

        # Make ElementTree from a string, then call suites_from_xml.
        xml_tree = ElementTree(XML(xml_str))
        spec_list = list(suites_from_xml(xml_tree, known_paths))

        self.assertEqual(len(names), len(directories),
                         msg="Internal test suite error: name and "+
                         "directories lists are different sizes!")

        self.assertEqual(len(spec_list), len(names),
                         msg="Wrong number of suite specs returned.")

        self.assertEqual(names,
                         [spec.name for spec in spec_list],
                         msg="Wrong suite name(s).")

        self.assertEqual(directories,
                         [spec.directories for spec in spec_list],
                         msg="Wrong suite path(s).")

        if labels is not None:
            self.assertEqual(labels,
                             [spec.labels for spec in spec_list],
                             msg="Wrong suite label(s).")

    def test_no_suites(self):
        """suites_from_xml output returns empty list for no matches."""

        xml_str = """
<root>
</root>
"""
        self.check_spec_list(xml_str, [], [])

    def test_single_suite(self):
        """suites_from_xml output is correct for a single match."""

        xml_str = """
<root>
 <suite name="suite1">
  <directory>/the/path</directory>
 </suite>
</root>
"""
        self.check_spec_list(xml_str, ["suite1"], [["/the/path"]])

    def test_multiple_suites(self):
        """suites_from_xml output is correct for multiple matches."""

        xml_str = """
<root>
 <suite name="suite1">
  <directory>/the/path</directory>
 </suite>
 <suite name="suite2">
  <directory>/other/path</directory>
 </suite>
</root>
"""
        self.check_spec_list(xml_str,
                             ["suite1", "suite2"], 
                             [["/the/path"], ["/other/path"]])

    def test_path_relative_to_known(self):
        """suites_from_xml handles a relative_to directory attribute."""
        from os.path import abspath

        xml_str = """
<root>
 <suite name="suite1">
  <directory relative_to="foo">path</directory>
 </suite>
</root>
"""
        self.check_spec_list(xml_str,
                             ["suite1"], 
                             [["/foodir/path"]],
                             known_paths={"foo": "/foodir"})

    def test_path_with_whitespace(self):
        """suites_from_xml handles a directory with whitespace added."""
        from os.path import abspath

        xml_str = """
<root>
 <suite name="suite1">
  <directory>
    /the/path
  </directory>
 </suite>
</root>
"""
        self.check_spec_list(xml_str, ["suite1"], [["/the/path"]])

    def test_path_with_label(self):
        """suites_from_xml handles a directory with a label correctly."""
        from os.path import abspath

        xml_str = """
<root>
 <suite name="suite1">
  <directory label="foo">/foo</directory>
 </suite>
</root>
"""
        self.check_spec_list(xml_str, ["suite1"], [["/foo"]],
                             labels=[["foo"]])


if __name__ == "__main__":
    unittest.main()
