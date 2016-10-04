"""Parser for an XML test list.

Public classes:
TestSuiteSpec - Specification of unit test suite directories.

Public functions:
suites_from_xml - Generator returning test suite descriptions.
"""
import os.path

__all__ = ("TestSuiteSpec", "suites_from_xml")

class TestSuiteSpec(object):

    """Specification for the location of a test suite.

    Public methods:
    __init__ - Constructor.
    __iter__ - Iterate over label/directory pairs.

    Public data:
    name - Name of the suite for logging output.
    labels - Labels for test suite directories.
    directories - Absolute paths to test suite directories.
    """

    UNLABELED_STRING = "UNLABELED"

    def __init__(self, name, labels, directories):
        """Constructor for TestSuiteSpec.

        Arguments:
        name - Name of the test suite.
        labels - Labels for each directory.
        directories - Path to the test suite.
        """

        assert (len(labels) == len(directories)), \
            "TestSuiteSpec: Number of spec labels and number of spec "+ \
            "directories do not match."

        self.name = name
        self.labels = []

        for label in labels:
            if label is not None:
                self.labels.append(label)
            else:
                self.labels.append(self.UNLABELED_STRING)

        self.directories = [os.path.abspath(directory)
                            for directory in directories]

    def __iter__(self):
        """Iterate over directories.

        Each iteration yields a (label, directory) pair.
        """
        return ( (l, d) for l, d in zip(self.labels, self.directories) )

def suites_from_xml(xml_tree, known_paths={}):
    """Generate test suite descriptions from XML.

    Returns a TestSuiteSpec for each suite description in the XML input.

    Arguments:
    xml_tree - A standard library ElementTree for an XML suite list.
    known_paths - A dict naming known locations on the system.
                  (Default is empty dict.)

    The expected layout of the xml tree is as follows:

    1. At the top of the tree are suite tags.
    2. Each suite tag has a name attribute.
    3. Each suite tag may have paths in directory subelements.
    4. A directory tag may have a "relative_to" attribute, naming a key in
       the known_paths dict.
    """

    elements = xml_tree.findall("suite")

    for elem in elements:
        labels = []
        directories = []
        for directory in elem.findall("directory"):
            path = directory.text.strip()
            if "relative_to" in directory.keys():
                relative_to_key = directory.get("relative_to")
                assert relative_to_key in known_paths, \
                    "suites_from_xml: Unrecognized relative_to attribute."
                path = os.path.join(known_paths[relative_to_key],
                                    path)
            directories.append(path)
            if "label" in directory.keys():
                labels.append(directory.get("label"))
            else:
                labels.append(None)
        yield TestSuiteSpec(elem.get("name"), labels, directories)
