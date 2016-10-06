"""Interface to `namelist_defaults.xml`.

This module contains only one class, `NamelistDefaults`, inheriting from
`GenericXML`.
"""

# Warnings we typically ignore.
# pylint:disable=invalid-name

# Disable warnings due to using `standard_module_setup`
# pylint:disable=wildcard-import,unused-wildcard-import

from CIME.XML.standard_module_setup import *
from CIME.XML.entry_id import EntryID

logger = logging.getLogger(__name__)

class NamelistDefaults(EntryID):

    """Class representing variable default values for a namelist.

    This class inherits from `GenericXML`, and supports most inherited methods;
    however, `get_resolved_value` is unsupported.

    Additional public methods:
    - add
    """

    def __init__(self, infile, attributes=None):
        """Construct a `NamelistDefaults` from an XML file."""
        super(NamelistDefaults, self).__init__(infile)
        self._attributes = attributes

    def add(self, infile):
        """Add the contents of an XML file to the defaults."""
        new_root = ET.parse(infile).getroot()
        for elem in new_root:
            self.root.append(elem)

    @staticmethod
    def _split_defaults_text(string):
        """Take a comma-separated list in a string, and split it into a list."""
        # Some trickiness here; we want to split items on commas, but not inside
        # quote-delimited strings. Stripping whitespace is also useful.
        value = []
        pos = 0
        delim = None
        for i, char in enumerate(string):
            if delim is None:
                # If not inside a string...
                if char in ('"', "'"):
                    # if we have a quote character, start a string.
                    delim = char
                elif char == ',':
                    # if we have a comma, this is a new value.
                    value.append(string[pos:i].strip())
                    pos = i+1
            else:
                # If inside a string, the only thing that can happen is the end
                # of the string.
                if char == delim:
                    delim = None
        value.append(string[pos:].strip())
        return value

    def get_value(self, item, attribute=None, resolved=False, subgroup=None):
        """Return the default value for the variable named `item`.

        The return value is a list of strings corresponding to the
        comma-separated list of entries for the value (length 1 for scalars). If
        there is no default value in the file, this returns `None`.
        """
        expect(not resolved, "This class does not support env resolution.")
        expect(subgroup is None, "This class does not support subgroups.")

        # #nodes = self.get_nodes(item.lower())
        # #node = self.get_node("entry", attribute={"id":item.lower()})
        # values = self.get_values(item.lower(), attribute=attribute, resolved=resolved, subgroup=subgroup)
        # print "DEBUG: values for item %s are %s" %(item,values)

        # Merge internal attributes with those passed in.
        all_attributes = {}
        if self._attributes is not None:
            all_attributes.update(self._attributes)
        if attribute is not None:
            all_attributes.update(attribute)
        value =  self.get_value_match(item.lower(), attributes=all_attributes)
        if value is not None:
            value =  self._split_defaults_text(value)
        return value

        # # Store nodes that match the attributes and their scores.
        # matches = []
        # for node in nodes:
        #     # For each node in the list start a score.
        #     score = 0
        #     for attribute in node.keys():
        #         # For each attribute, add to the score.
        #         score += 1
        #         # If some attribute is specified that we don't know about,
        #         # or the values don't match, it's not a match we want.
        #         if attribute not in all_attributes or \
        #            all_attributes[attribute] != node.get(attribute):
        #             score = -1
        #             break

        #     # Add valid matches to the list.
        #     if score >= 0:
        #         matches.append((score, node))

        # if not matches:
        #     return None

        # # Get maximum score using custom `key` function, extract the node.
        # _, node = max(matches, key=lambda x: x[0])
        # if node.text is None:
        #     return ['']
        # return self._split_defaults_text(node.text)

    # While there are environment variable references in the file at times, they
    # usually involve env files that we don't know about at this low level. So
    # we punt on this issue and make the caller resolve these, if necessary.
    def get_resolved_value(self, raw_value):
        """This function is not implemented."""
        raise TypeError, \
            "NamelistDefaults does not support `get_resolved_value`."
