"""Interface to `namelist_definition.xml`.

This module contains only one class, `NamelistDefinition`, inheriting from
`EntryID`.
"""

# Warnings we typically ignore.
# pylint:disable=invalid-name

# Disable warnings due to using `standard_module_setup`
# pylint:disable=wildcard-import,unused-wildcard-import

import re
import collections

from CIME.namelist import fortran_namelist_base_value, \
    is_valid_fortran_namelist_literal, character_literal_to_string, \
    expand_literal_list, Namelist, get_fortran_name_only

from CIME.XML.standard_module_setup import *
from CIME.XML.entry_id import EntryID
from CIME.XML.files import Files

logger = logging.getLogger(__name__)

_array_size_re = re.compile(r'^(?P<type>[^(]+)\((?P<size>[^)]+)\)$')

class CaseInsensitiveDict(dict):

    """Basic case insensitive dict with strings only keys.
        From https://stackoverflow.com/a/27890005 """

    proxy = {}

    def __init__(self, data):
        dict.__init__(self)
        self.proxy = dict((k.lower(), k) for k in data)
        for k in data:
            self[k] = data[k]

    def __contains__(self, k):
        return k.lower() in self.proxy

    def __delitem__(self, k):
        key = self.proxy[k.lower()]
        super(CaseInsensitiveDict, self).__delitem__(key)
        del self.proxy[k.lower()]

    def __getitem__(self, k):
        key = self.proxy[k.lower()]
        return super(CaseInsensitiveDict, self).__getitem__(key)

    def get(self, k, default=None):
        return self[k] if k in self else default

    def __setitem__(self, k, v):
        super(CaseInsensitiveDict, self).__setitem__(k, v)
        self.proxy[k.lower()] = k

class NamelistDefinition(EntryID):

    """Class representing variable definitions for a namelist.
    This class inherits from `EntryID`, and supports most inherited methods;
    however, `set_value` is unsupported.

    Additional public methods:
    - dict_to_namelist.
    - is_valid_value
    - validate
    """

    def __init__(self, infile, files=None):
        """Construct a `NamelistDefinition` from an XML file."""

        # if the file is invalid we may not be able to check the version
        # but we need to do it this way until we remove the version 1 files
        schema = None
        if files is None:
            files = Files()
        schema = files.get_schema("NAMELIST_DEFINITION_FILE")
        expect(os.path.isfile(infile), "File {} does not exist".format(infile))
        super(NamelistDefinition, self).__init__(infile, schema=schema)

        self._attributes = {}
        self._entry_nodes = []
        self._entry_ids = []
        self._valid_values = {}
        self._entry_types = {}
        self._group_names = CaseInsensitiveDict({})
        self._nodes = {}

    def set_nodes(self, skip_groups=None):
        """
        populates the object data types for all nodes that are not part of the skip_groups array
        returns nodes that do not have attributes of `skip_default_entry` or `per_stream_entry`
        """
        default_nodes = []
        for node in self.get_children("entry"):
            name = self.get(node, "id")
            skip_default_entry = self.get(node, "skip_default_entry") == "true"
            per_stream_entry = self.get(node, "per_stream_entry") == "true"
            set_node_values = False
            if skip_groups:
                group_name = self._get_group_name(node)
                if not group_name in skip_groups:
                    self._entry_nodes.append(node)
                    set_node_values = True
                    if not skip_default_entry and not per_stream_entry:
                        default_nodes.append(node)
            else:
                self._entry_nodes.append(node)
                set_node_values = True
                if not skip_default_entry and not per_stream_entry:
                    default_nodes.append(node)
            if set_node_values:
                self._entry_nodes.append(node)
                self._entry_ids.append(name)
                self._nodes[name] = node
                self._entry_types[name] = self._get_type(node)
                self._valid_values[name] = self._get_valid_values(node)
                self._group_names[name] = self._get_group_name(node)
        return default_nodes

    def _get_group_name(self, node=None):
        if self.get_version() == 1.0:
            group = self.get(node, 'group')
        elif self.get_version() >= 2.0:
            group = self.get_element_text("group", root=node)
        return(group)

    def _get_type(self, node):
        if self.get_version() == 1.0:
            type_info = self.get(node, 'type')
        elif self.get_version() >= 2.0:
            type_info = self._get_type_info(node)
        return(type_info)

    def _get_valid_values(self, node):
        # The "valid_values" attribute is not required, and an empty string has
        # the same effect as not specifying it.
        # Returns a list from a comma seperated string in xml
        valid_values = ''
        if self.get_version() == 1.0:
            valid_values = self.get(node, 'valid_values')
        elif self.get_version() >= 2.0:
            valid_values = self._get_node_element_info(node, "valid_values")
        if valid_values == '':
            valid_values = None
        if valid_values is not None:
            valid_values = valid_values.split(',')
        return valid_values

    def get_group(self, name):
        return self._group_names[name]

    def add_attributes(self, attributes):
        self._attributes = attributes

    def get_entry_nodes(self):
        return self._entry_nodes

    def get_per_stream_entries(self):
        entries = []
        nodes = self.get_children("entry")
        for node in nodes:
            per_stream_entry = self.get(node, "per_stream_entry") == "true"
            if per_stream_entry:
                entries.append(self.get(node, "id"))
        return entries

    # Currently we don't use this object to construct new files, and it's no
    # good for that purpose anyway, so stop this function from being called.
    def set_value(self, vid, value, subgroup=None, ignore_type=True):
        """This function is not implemented."""
        raise TypeError("NamelistDefinition does not support `set_value`.")

    def get_value_match(self, vid, attributes=None, exact_match=True, entry_node=None):
        """Return the default value for the variable named `vid`.

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

        if entry_node is None:
            entry_node = self._nodes[vid]
        value = super(NamelistDefinition, self).get_value_match(vid.lower(),attributes=all_attributes, exact_match=exact_match,
                                                                entry_node=entry_node)
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

    def split_type_string(self, name):
        """Split a 'type' attribute string into its component parts.

        The `name` argument is the variable name.
        This is used for error reporting purposes.

        The return value is a tuple consisting of the type itself, a length
        (which is an integer for character variables, otherwise `None`), and the
        size of the array (which is 1 for scalar variables).
        """
        type_string = self._entry_types[name]

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
                       "In namelist definition, variable {} had the non-integer string {!r} specified as an array size.".format(name, size_string))
        else:
            size = 1

        # Separate into a type and an optional length.
        type_, star, length = type_string.partition('*')
        if star == '*':
            # Length allowed only for character variables.
            expect(type_ == 'character',
                   "In namelist definition, length specified for non-character "
                   "variable {}.".format(name))
            # Check that the length is actually an integer, to make the error
            # message a bit cleaner if the xml input is bad.
            try:
                max_len = int(length)
            except ValueError:
                expect(False,
                       "In namelist definition, character variable {} had the non-integer string {!r} specified as a length.".format(name, length))
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
        # Separate into a type, optional length, and optional size.
        type_, max_len, size = self.split_type_string(name)
        invalid = []

        # Check value against type.
        for scalar in value:
            if not is_valid_fortran_namelist_literal(type_, scalar):
                invalid.append(scalar)
        if len(invalid) > 0:
            logger.warning("Invalid values {}".format(invalid))
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
        valid_values = self._valid_values[name]
        if valid_values is not None:
            expect(type_ in ('integer', 'character'),
                   "Found valid_values attribute for variable {} with type {}, but valid_values only allowed for character and integer variables.".format(name, type_))
            if type_ == 'integer':
                compare_list = [int(vv) for vv in valid_values]
            else:
                compare_list = valid_values
            for scalar in canonical_value:
                if scalar not in compare_list:
                    invalid.append(scalar)
            if len(invalid) > 0:
                logger.warning("Invalid values {}".format(invalid))
                return False

        # Check size of input array.
        if len(expand_literal_list(value)) > size:
            expect(False, "Value index exceeds variable size for variable {}, allowed array length is {} value array size is {}".format(name, size, len(expand_literal_list(value))))
        return True

    def _expect_variable_in_definition(self, name, variable_template):
        """Used to get a better error message for an unexpected variable.
             case insensitve match"""

        expect(name in self._entry_ids,
               (variable_template + " is not in the namelist definition.").format(str(name)))

    def _user_modifiable_in_variable_definition(self, name):
        # Is name user modifiable?
        node = self.get_optional_child("entry", attributes={'id': name})
        user_modifiable_only_by_xml = self.get(node, 'modify_via_xml')
        if user_modifiable_only_by_xml is not None:
            expect(False,
                   "Cannot change {} in user_nl file: set via xml variable {}".format(name, user_modifiable_only_by_xml))
        user_cannot_modify = self.get(node, 'cannot_modify_by_user_nl')
        if user_cannot_modify is not None:
            expect(False,
                   "Cannot change {} in user_nl file: {}".format(name, user_cannot_modify))
    def _generate_variable_template(self, filename):
        # Improve error reporting when a file name is provided.
        if filename is None:
            variable_template = "Variable {!r}"
        else:
            # for the next step we want the name of the original user_nl file not the internal one
            # We do this by extracting the component name from the filepath string
            if "Buildconf" in filename and "namelist_infile" in filename:
                msgfn = "user_nl_" + (filename.split(os.sep)[-2])[:-4]
            else:
                msgfn = filename
            variable_template = "Variable {!r} from file " + repr(str(msgfn))
        return variable_template

    def validate(self, namelist,filename=None):
        """Validate a namelist object against this definition.

        The optional `filename` argument can be used to assist in error
        reporting when the namelist comes from a specific, known file.
        """
        variable_template = self._generate_variable_template(filename)

        # Iterate through variables.
        for group_name in namelist.get_group_names():
            for variable_name in namelist.get_variable_names(group_name):
                # Check that the variable is defined...
                qualified_variable_name = get_fortran_name_only(variable_name)
                self._expect_variable_in_definition(qualified_variable_name, variable_template)

                # Check if can actually change this variable via filename change
                if filename is not None:
                    self._user_modifiable_in_variable_definition(qualified_variable_name)

                # and has the right group name...
                var_group = self.get_group(qualified_variable_name)
                expect(var_group == group_name,
                       (variable_template + " is in a group named {!r}, but should be in {!r}.").format(str(variable_name), str(group_name), str(var_group)))

                # and has a valid value.
                value = namelist.get_variable_value(group_name, variable_name)
                expect(self.is_valid_value(qualified_variable_name, value),
                       (variable_template + " has invalid value {!r}.").format(str(variable_name), [str(scalar) for scalar in value]))

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
        variable_template = self._generate_variable_template(filename)
        groups = {}
        for variable_name in dict_:
            variable_lc = variable_name.lower()
            qualified_varname = get_fortran_name_only(variable_lc)
            self._expect_variable_in_definition(qualified_varname, variable_template)
            group_name = self.get_group(qualified_varname)
            expect (group_name is not None, "No group found for var {}".format(variable_lc))
            if group_name not in groups:
                groups[group_name] = collections.OrderedDict()
            groups[group_name][variable_lc] = dict_[variable_name]
        return Namelist(groups)

    def get_input_pathname(self, name):
        node = self._nodes[name]
        if self.get_version() == 1.0:
            input_pathname = self.get(node, 'input_pathname')
        elif self.get_version() >= 2.0:
            input_pathname = self._get_node_element_info(node, "input_pathname")
        return(input_pathname)

    # pylint: disable=arguments-differ
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
