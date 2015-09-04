#!/usr/bin/env python
"""Unit tests for the xml_utils module.

Public classes:
TestBestMatch - best_match tests.
TestAllMatches - all_matches tests.
TestElementsToDict - elements_to_dict tests.
"""

import unittest
from xml.etree.ElementTree import XML, ElementTree

from xml_utils import best_match, all_matches, elements_to_dict

__all__ = ("TestBestMatch", "TestAllMatches", "TestElementsToDict")

class TestBestMatch(unittest.TestCase):

    """Tests for the best_match function."""

    def setUp(self):
        """Create an ElementTree object as test data.

        This is a complex set of data, but the tests on it are
        comparatively simple. Some of the data (such as the empty "foo"
        tags) are present purely to prove that they are never matched.
        """

        root = XML("""
<root>
 <foo a="1">
  <data2>a</data2>
  <data3>a</data3>
 </foo>
 <foo>
  <data1>generic</data1>
  <data2>generic</data2>
  <data3>generic</data3>
 </foo>
 <foo b="2">
  <data2>b</data2>
  <data3>b</data3>
 </foo>
 <foo a="1" b="2">
  <data3>ab</data3>
 </foo>
 <foo a="2" />
 <foo b="1" />
 <foo a="1" b="1" />
 <foo a="2" b="1" />
 <bar a="2">
  <data1>abar</data1>
 </bar>
 <bar b="1">
  <data1>bbar</data1>
 </bar>
</root>
""")

        self.xml_tree = ElementTree(root)

    def test_no_match(self):
        """best_match returns None when no paths match."""

        self.assertTrue(best_match(self.xml_tree, "invalid")
                        is None)

    def test_simple_match(self):
        """best_match can find the only match, when no attributes are used."""

        self.assertEqual(
            "generic",
            best_match(self.xml_tree, "foo/data1").text
            )

    def test_no_match_attr(self):
        """best_match returns None when there are no attribute matches."""

        self.assertTrue(best_match(self.xml_tree, "bar/data1")
                        is None)
        self.assertTrue(best_match(self.xml_tree, "bar/data1", {"a": "1"})
                        is None)

    def test_match_with_attr(self):
        """best_match returns the only path match, with matching attribute."""

        self.assertEqual(
            "abar",
            best_match(self.xml_tree, "bar/data1", {"a": "2"}).text
            )

    def test_match_on_attr(self):
        """best_match returns the only path and attribute match."""

        self.assertEqual(
            "b",
            best_match(self.xml_tree, "foo/data2",
                       {"a": "2", "b": "2"}).text
            )

    def test_match_most_specific(self):
        """best_match returns the most specific match for each path."""

        self.assertEqual(
            "generic",
            best_match(self.xml_tree, "foo/data1",
                       {"a": "1", "b": "2"}).text
            )
        self.assertEqual(
            "ab",
            best_match(self.xml_tree, "foo/data3",
                       {"a": "1", "b": "2"}).text
            )

    def test_match_first(self):
        """best_match returns the first matching entry."""

        self.assertEqual(
            "abar",
            best_match(self.xml_tree, "bar/data1",
                       {"a": "2", "b": "1"}).text
            )


class TestAllMatches(unittest.TestCase):

    """Tests for the all_matches function."""

    def setUp(self):
        """Create an ElementTree for use as test data."""

        root = XML("""
<root>
 <test1>
  <data>the_text</data>
 </test1>
 <test2>
  <data>first_text</data>
  <data>second_text</data>
 </test2>
 <test3>
 </test3>
 <test4>
  <data invalid="true">first_text</data>
  <data valid="false">second_text</data>
  <data valid="true" invalid="true">second_text</data>
 </test4>
 <test5>
  <data invalid="false">first_text</data>
  <data valid="true">second_text</data>
  <data>third_text</data>
  <data valid="true" invalid="false">fourth_text</data>
  <data valid="true" invalid="true">bad_text</data>
 </test5>
 <test6>
  <data ignore="ignore this">the_text</data>
 </test6>
</root>
""")

        self.xml_tree = ElementTree(root)

    def test_one_match(self):
        """all_matches returns one element successfully."""
        self.assertEqual(
            [e.text for e in all_matches(self.xml_tree, "test1/data")],
            ["the_text"],
            )

    def test_multiple_matches(self):
        """all_matches returns many elements successfully.

        Elements must be returned in the right order as well.
        """
        self.assertEqual(
            [e.text for e in all_matches(self.xml_tree, "test2/data")],
            ["first_text", "second_text"],
            )

    def test_no_matches(self):
        """all_matches returns no elements if none are found."""
        self.assertEqual(
            [e.text for e in all_matches(self.xml_tree, "test3/data")],
            [],
            )

    def test_no_attr_matches(self):
        """all_matches returns no elements if none match attributes."""
        self.assertEqual(
            [e.text for e in all_matches(self.xml_tree, "test4/data",
                                         {"valid": "true"})],
            [],
            )

    def test_attr_matches(self):
        """all_matches returns elements that match attributes."""
        self.assertEqual(
            [e.text for e in
             all_matches(self.xml_tree, "test5/data",
                         {"valid": "true", "invalid": "false"})],
            ["first_text", "second_text", "third_text", "fourth_text"],
            )

    def test_ignore_attribute(self):
        """all_matches can ignore "extra" attributes."""
        self.assertEqual(
            [e.text for e in
             all_matches(self.xml_tree, "test6/data",
                         ignore=["ignore"])],
            ["the_text"]
            )


# In true TDD form, these tests are far longer than the code they apply to!
# Of course, when they were written, the first version of the code was
# longer.
class TestElementsToDict(unittest.TestCase):

    """Tests for the elements_to_dict function."""

    @staticmethod
    def string_to_elements(xml_str, ignore="key"):
        """Converts an XML string to a list of its "ENV" elements."""
        return all_matches(ElementTree(XML(xml_str)), "ENV",
                           ignore=ignore)

    def test_no_elements(self):
        """elements_to_dict with no elements produces an empty dict."""
        self.assertEqual(elements_to_dict([]), {})

    def test_one_element(self):
        """elements_to_dict can handle one element."""
        elements = self.string_to_elements("""
<root>
 <ENV key="foo">bar</ENV>
</root>
""")
        self.assertEqual(elements_to_dict(elements), {"foo": "bar"})

    def test_multiple_elements(self):
        """elements_to_dict can handle multiple elements."""
        elements = self.string_to_elements("""
<root>
 <ENV key="foo1">bar1</ENV>
 <ENV key="foo2">bar2</ENV>
</root>
""")
        self.assertEqual(elements_to_dict(elements),
                         {"foo1": "bar1", "foo2": "bar2"})

    def test_key_attr(self):
        """elements_to_dict uses key_attr as key attribute."""
        key_attr = "name"
        xml_str = """
<root>
 <ENV name="foo">bar</ENV>
</root>
"""
        elements = self.string_to_elements(xml_str, ignore=key_attr)
        self.assertEqual(elements_to_dict(elements, key_attr=key_attr),
                         {"foo": "bar"})

    def test_multiple_key(self):
        """elements_to_dict lets later values overwrite earlier ones.

        This might not actually be the best behavior, but it's the easiest
        to implement; this test is really just to prevent *accidentally*
        changing the design.
        """
        elements = self.string_to_elements("""
<root>
 <ENV key="foo">bar1</ENV>
 <ENV key="foo">bar2</ENV>
</root>
""")
        self.assertEqual(elements_to_dict(elements),
                         {"foo": "bar2"})

    def test_keyless_match(self):
        """elements_to_dict with a match with no key returns empty dict.

        It might be better to raise an exception; this is a robustness and
        simplicity vs. correctness tradeoff.
        """
        elements = self.string_to_elements("""
<root>
 <ENV>bar</ENV>
</root>
""")
        self.assertEqual(elements_to_dict(elements), {})


if __name__ == "__main__":
    unittest.main()
