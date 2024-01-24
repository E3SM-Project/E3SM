#!/usr/bin/env python3

import os
import unittest
import urllib.request
import xml.etree.ElementTree as ET


class testNamelistDefaultsScream(unittest.TestCase):
    def setUp(self):
        """
        Set up the environment for the test by setting the DIN_LOC_ROOT
        environment variable. Parse the 'namelist_defaults_scream.xml' 
        file and extract the files of interest based on the DIN_LOC_ROOT
        variable. Assign the extracted files to the 'my_files' attribute
        of the test instance.
        """
        os.environ["DIN_LOC_ROOT"] = "https://web.lcrc.anl.gov/public/e3sm/inputdata/"

        tree = ET.parse('../namelist_defaults_scream.xml')
        root = tree.getroot()

        file_of_interest = [
            child.text for child in root.findall(".//")
            if child.text and child.text.startswith("${DIN_LOC_ROOT}")
        ]

        self.my_files = [file.replace("${DIN_LOC_ROOT}/", "")
                         for file in file_of_interest]

    def test_retrieve_files(self):
        """
        Function to test the retrieval of files. 
        """
        for i_file in range(len(self.my_files)):
            with self.subTest(i=i_file):
                request_return = urllib.request.urlretrieve(
                    f"{os.environ['DIN_LOC_ROOT']}{self.my_files[i_file]}"
                )
                self.assertIsInstance(request_return, tuple)

    def test_expected_fail(self):
        """
        A test case to verify the expected failure when encountering
        an HTTP error using urllib.
        """
        with self.assertRaises(urllib.error.HTTPError):
            some_phony_file = f"{self.my_files[5][:-5]}some_phony_file.nc"
            urllib.request.urlretrieve(
                f"{os.environ['DIN_LOC_ROOT']}{some_phony_file}"
            )


if __name__ == '__main__':
    unittest.main()
