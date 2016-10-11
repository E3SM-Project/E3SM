"""Interface to `namelist_definition.xml`.

This module contains only one class, `NamelistDefinition`, inheriting from
`GenericXML`.
"""

# Warnings we typically ignore.
# pylint:disable=invalid-name

# Disable warnings due to using `standard_module_setup`
# pylint:disable=wildcard-import,unused-wildcard-import

import re

from CIME.namelist import fortran_namelist_base_value, \
    is_valid_fortran_namelist_literal, character_literal_to_string, \
    expand_literal_list, Namelist

from CIME.XML.standard_module_setup import *
from CIME.XML.entry_id import EntryID
from CIME.utils import get_cime_root

logger = logging.getLogger(__name__)

_array_size_re = re.compile(r'^(?P<type>[^(]+)\((?P<size>[^)]+)\)$')

class NamelistDefinition(EntryID):

    """Class representing variable definitions for a namelist.
    This class inherits from `EntryID`, and supports most inherited methods;
    however, `set_value` is unsupported.

    Additional public methods:
    - dict_to_namelist.
    - is_valid_value
    - validate
    """

    def __init__(self, infile, attributes=None):
        """Construct a `NamelistDefinition` from an XML file."""
        super(NamelistDefinition, self).__init__(infile)
        self._attributes = attributes
        # if the file is invalid we may not be able to check the version
        # but we need to do it this way until we remove the version 1 files
#        if self._get_version() == "2.0":
#            cimeroot = get_cime_root()
#            schema = os.path.join(cimeroot,"cime_config","xml_schemas","entry_id_namelist.xsd")
#            self.validate_xml_file(infile, schema)


    def _get_version(self):
        version = self.root.get("version")
        return version

    def get_entries(self):
        """Return all variables in the namelist definition file
        that do not have attributes of skip_default_entry or per_stream_entry
        """
        entries = []
        nodes = self.get_nodes("entry")
        for node in nodes:
            skip_default_entry = node.get("skip_default_entry")
            per_stream_entry = node.get("per_stream_entry")
            if not skip_default_entry and not per_stream_entry:
                entries.append(node.get("id"))
        return entries

    def get_per_stream_entries(self):
        entries = []
        nodes = self.get_nodes("entry")
        for node in nodes:
            per_stream_entry = node.get("per_stream_entry")
            if per_stream_entry:
                entries.append(node.get("id"))
        return entries

    # Currently we don't use this object to construct new files, and it's no
    # good for that purpose anyway, so stop this function from being called.
    def set_value(self, vid, value, subgroup=None, ignore_type=True):
        """This function is not implemented."""
        raise TypeError, \
            "NamelistDefinition does not support `set_value`."

    def get_valid_values(self, name):
        # The "valid_values" attribute is not required, and an empty string has
        # the same effect as not specifying it.
        # Returns a list from a comma seperated string in xml
        valid_values = ''
        elem = self.get_optional_node("entry", attributes={'id': name})
        if self._get_version() == "1.0":
            valid_values = elem.get('valid_values')
        elif self._get_version() == "2.0":
            valid_values = self._get_node_element_info(elem, "valid_values")
        if valid_values == '':
            valid_values = None
        if valid_values is not None:
            valid_values = valid_values.split(',')
        return valid_values

    def get_value_match(self, item, attributes=None, exact_match=True):
        """Return the default value for the variable named `item`.

        The return value is a list of strings corresponding to the
        comma-separated list of entries for the value (length 1 for scalars). If
        there is no default value in the file, this returns `None`.
        """
        # Merge internal attributes with those passed in.
        all_attributes = {}
        if self._attributes is not None:
            all_attributes.update(self._attributes)
        if attributes is not None:
            all_attributes.update(attributes)

        value = super(NamelistDefinition, self).get_value_match(item.lower(),attributes=all_attributes, exact_match=exact_match)

        if value is None:
            value = ''
        else:
            value =  self._split_defaults_text(value)

        return value

    @staticmethod
    def _split_defaults_text(string):
        """Take a comma-separated list in a string, and split it into a list."""
        # Some trickiness here; we want to split items on commas, but not inside
        # quote-delimited strings. Stripping whitespace is also useful.
        value = []
        if len(string):
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

    @staticmethod
    def split_type_string(name, type_string):
        """Split a 'type' attribute string into its component parts.

        The `name` argument is the variable name associated with this type
        string. It is used for error reporting purposes.

        The return value is a tuple consisting of the type itself, a length
        (which is an integer for character variables, otherwise `None`), and the
        size of the array (which is 1 for scalar variables).

        This method also checks to ensure that the input `type_string` is valid.
        """
        # 'char' is frequently used as an abbreviation of 'character'.
        type_string = type_string.replace('char', 'character')
        # Separate into a size and the rest of the type.
        size_match = _array_size_re.search(type_string)
        if size_match:
            type_string = size_match.group('type')
            size_string = size_match.group('size')
            try:
                size = int(size_string)
            except ValueError:
                expect(False,
                       "In namelist definition, variable %s had the "
                       "non-integer string %r specified as an array size." %
                       (name, size_string))
        else:
            size = 1
        # Separate into a type and an optional length.
        type_, star, length = type_string.partition('*')
        if star == '*':
            # Length allowed only for character variables.
            expect(type_ == 'character',
                   "In namelist definition, length specified for non-character "
                   "variable %s." % name)
            # Check that the length is actually an integer, to make the error
            # message a bit cleaner if the xml input is bad.
            try:
                max_len = int(length)
            except ValueError:
                expect(False,
                       "In namelist definition, character variable %s had the "
                       "non-integer string %r specified as a length." %
                       (name, length))
        else:
            max_len = None
        return type_, max_len, size

    @staticmethod
    def _canonicalize_value(type_, value):
        """Create 'canonical' version of a value for comparison purposes."""
        canonical_value = [fortran_namelist_base_value(scalar)
                           for scalar in value]
        canonical_value = [scalar for scalar in canonical_value if scalar != '']
        if type_ == 'character':
            canonical_value = [character_literal_to_string(scalar)
                               for scalar in canonical_value]
        elif type_ == 'integer':
            canonical_value = [int(scalar) for scalar in canonical_value]
        return canonical_value

    def is_valid_value(self, name, value):
        """Determine whether a value is valid for the named variable.

        The `value` argument must be a list of strings formatted as they would
        appear in the namelist (even for scalar variables, in which case the
        length of the list is always 1).
        """
        name = name.lower()
        # Separate into a type, optional length, and optional size.
        type_, max_len, size = self.split_type_string(name, self.get_type_info(name))

        # Check value against type.
        for scalar in value:
            if not is_valid_fortran_namelist_literal(type_, scalar):
                return False

        # Now that we know that the strings as input are valid Fortran, do some
        # canonicalization for further checks.
        canonical_value = self._canonicalize_value(type_, value)

        # Check maximum length (if applicable).
        if max_len is not None:
            for scalar in canonical_value:
                if len(scalar) > max_len:
                    return False

        # Check valid value constraints (if applicable).
        valid_values = self.get_valid_values(name)
        if valid_values is not None:
            expect(type_ in ('integer', 'character'),
                   "Found valid_values attribute for variable %s with "
                   "type %s, but valid_values only allowed for character "
                   "and integer variables." % (name, type_))
            if type_ == 'integer':
                compare_list = [int(vv)
                                for vv in valid_values]
            else:
                compare_list = valid_values
            for scalar in canonical_value:
                if scalar not in compare_list:
                    return False

        # Check size of input array.
        if len(expand_literal_list(value)) > size:
            return False
        return True

    def _expect_variable_in_definition(self, name, variable_template):
        """Used to get a better error message for an unexpected variable."""
        node = self.get_optional_node("entry", attributes={'id': name})
        expect(node is not None,
               (variable_template + " is not in the namelist definition.") %
               str(name))

    def _user_modifiable_in_variable_definition(self, name):
        # Is name user modifiable?
        node = self.get_optional_node("entry", attributes={'id': name})
        user_modifiable = node.get('modify_via_xml')
        if user_modifiable is not None:
            expect(False,
                   "Cannot change %s in user_nl_xxx file, %s" %(name, user_modifiable))

    def validate(self, namelist, filename=None):
        """Validate a namelist object against this definition.

        The optional `filename` argument can be used to assist in error
        reporting when the namelist comes from a specific, known file.
        """
        # Improve error reporting when a file name is provided.
        if filename is None:
            variable_template = "Variable %r"
        else:
            variable_template = "Variable %r from file " + repr(str(filename))

        # Iterate through variables.
        for group_name in namelist.get_group_names():
            for variable_name in namelist.get_variable_names(group_name):
                # Check that the variable is defined...
                self._expect_variable_in_definition(variable_name, variable_template)

                # Check if can actually change this variable via filename change
                if filename is not None:
                    self._user_modifiable_in_variable_definition(variable_name)

                # and has the right group name...
                var_group = self.get_group_name(variable_name)
                expect(var_group == group_name,
                       (variable_template + " is in a group named %r, but "
                        "should be in %r.") %
                       (str(variable_name), str(group_name),
                        str(var_group)))

                # and has a valid value.
                value = namelist.get_variable_value(group_name, variable_name)
                expect(self.is_valid_value(variable_name, value),
                       (variable_template + " has invalid value %r.") %
                       (str(variable_name), [str(scalar) for scalar in value]))

    def dict_to_namelist(self, dict_, filename=None):
        """Converts a dictionary of name-value pairs to a `Namelist`.

        The input is assumed to be similar to the output of `parse` when
        `groupless=True` is set. This function uses the namelist definition file
        to look up the namelist group associated with each variable, and uses
        this information to create a true `Namelist` object.

        The optional `filename` argument can be used to assist in error
        reporting when the namelist comes from a specific, known file.
        """
        # Improve error reporting when a file name is provided.
        if filename is None:
            variable_template = "Variable %rs"
        else:
            variable_template = "Variable %r from file " + repr(str(filename))
        groups = {}
        for variable_name in dict_:
            variable_lc = variable_name.lower()
            self._expect_variable_in_definition(variable_lc, variable_template)
            group_name = self.get_group_name(variable_lc)
            expect (group_name is not None, "No group found for var %s"%variable_lc)
            if group_name not in groups:
                groups[group_name] = {}
            groups[group_name][variable_lc] = dict_[variable_name]
        return Namelist(groups)

    def get_input_pathname(self, name):
        elem = self.get_optional_node("entry", attributes={'id': name})

        if self._get_version() == "1.0":
            input_pathname = elem.get('input_pathname')
        elif self._get_version() == "2.0":
            input_pathname = self._get_node_element_info(elem, "input_pathname")
        return(input_pathname)

    def get_type_info(self, name):
        elem = self.get_optional_node("entry", attributes={'id': name})

        if self._get_version() == "1.0":
            type_info = elem.get('type')
        elif self._get_version() == "2.0":
            type_info = self._get_type_info(elem)
        return(type_info)

    def get_group_name(self, name):
        elem = self.get_optional_node("entry", attributes={'id': name})

        if self._get_version() == "1.0":
            group = elem.get('group')
        elif self._get_version() == "2.0":
            group = self._get_node_element_info(elem, "group")
        return(group)

    def get_default_value(self, item, attribute=None):
        """Return the default value for the variable named `item`.

        The return value is a list of strings corresponding to the
        comma-separated list of entries for the value (length 1 for scalars). If
        there is no default value in the file, this returns `None`.
        """
        # Merge internal attributes with those passed in.
        all_attributes = {}
        if self._attributes is not None:
            all_attributes.update(self._attributes)
        if attribute is not None:
            all_attributes.update(attribute)

        value = self.get_value_match(item.lower(), all_attributes, True)
        return self._split_defaults_text(value)
