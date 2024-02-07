#!/usr/bin/env python3

import os
import http
import pathlib
import unittest
import urllib.request
import xml.etree.ElementTree as ET


class testNamelistDefaultsScream(unittest.TestCase):
    def setUp(self):
        """
        Set up the environment for the test by setting the DIN_LOC_ROOT
        environment variable. Parse the 'namelist_defaults_scream.xml' 
        file and extract the files of interest based on the DIN_LOC_ROOT
        variable or the array(file) type. Assign the extracted files
        to the 'my_files' attribute of the test instance.
        """

        os.environ["DIN_LOC_ROOT"] = "https://web.lcrc.anl.gov/public/e3sm/inputdata/"

        scream_defaults_path = pathlib.Path(__file__)
        tree = ET.parse(f"{scream_defaults_path.parent.parent}/namelist_defaults_scream.xml")
        root = tree.getroot()

        files_of_interest = [
            child.text for child in root.findall(".//")
            if child.text and child.text.startswith("${DIN_LOC_ROOT}")
        ]

        more_files_of_interest = [
            child.text for child in root.findall(".//")
            if child.text and "type" in child.attrib.keys() and child.attrib["type"]=="array(file)"
        ]

        files_of_interest.extend(
            text.strip() for text_list in more_files_of_interest for text in text_list.split(",")
            if text.strip().startswith("${DIN_LOC_ROOT}")
        )

        self.my_files = [
            file.replace("${DIN_LOC_ROOT}/", "")
            for file in files_of_interest
        ]

        self.my_lines = []
        with open(
            f"{scream_defaults_path.parent.parent}/namelist_defaults_scream.xml",
            "r"
        ) as the_file:
            for a_line in the_file:
                self.my_lines.append(a_line)

    def test_ascii_lines(self):
        """
        Test that all lines are ASCII
        """

        for i_line, a_line in enumerate(self.my_lines):
            with self.subTest(i_line=i_line):
                self.assertTrue(
                    a_line.isascii(),
                    msg=f"\nERROR! This line is not ASCII!\n{a_line}"
                )

    def test_opening_files(self):
        """
        Test the opening of files from the inputdata server.
        """

        for i_file in range(len(self.my_files)):
            with self.subTest(i_file=i_file):
                try:
                    request_return = urllib.request.urlopen(
                        f"{os.environ['DIN_LOC_ROOT']}{self.my_files[i_file]}"
                    )
                    self.assertIsInstance(request_return, http.client.HTTPResponse)
                except urllib.error.HTTPError:
                    file_name = f"{os.environ['DIN_LOC_ROOT']}{self.my_files[i_file]}"
                    self.assertTrue(
                        False,
                        msg=f"\nERROR! This file doesn't exist!\n{file_name}"
                    )

    def test_expected_fail(self):
        """
        Test an expected failure by manipulating the file name.
        """

        with self.assertRaises(urllib.error.HTTPError):
            some_phony_file = f"{self.my_files[5][:-5]}some_phony_file.nc"
            urllib.request.urlopen(
                f"{os.environ['DIN_LOC_ROOT']}{some_phony_file}"
            )


if __name__ == '__main__':
    unittest.main()
