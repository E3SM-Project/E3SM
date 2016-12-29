"""Class for generating component namelists."""

# Typically ignore this.
# pylint: disable=invalid-name

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import

import datetime
import re

from CIME.XML.standard_module_setup import *
from CIME.namelist import Namelist, parse, \
    character_literal_to_string, string_to_character_literal, \
    expand_literal_list, compress_literal_list, merge_literal_lists
from CIME.XML.namelist_definition import NamelistDefinition
from CIME.utils import expect

logger = logging.getLogger(__name__)

_var_ref_re = re.compile(r"\$(\{)?(?P<name>\w+)(?(1)\})")

_ymd_re = re.compile(r"%(?P<digits>[1-9][0-9]*)?y(?P<month>m(?P<day>d)?)?")

_stream_file_template = """
<dataSource>
   GENERIC
</dataSource>
<domainInfo>
  <variableNames>
     {domain_varnames}
  </variableNames>
  <filePath>
     {domain_filepath}
  </filePath>
  <fileNames>
     {domain_filenames}
  </fileNames>
</domainInfo>
<fieldInfo>
   <variableNames>
     {data_varnames}
   </variableNames>
   <filePath>
     {data_filepath}
   </filePath>
   <fileNames>
    {data_filenames}
   </fileNames>
   <offset>
      {offset}
   </offset>
</fieldInfo>
"""

class NamelistGenerator(object):

    """Utility class for generating namelists for a given component."""

    _streams_variables = []

    def __init__(self, case, infiles, #pylint:disable=too-many-arguments
                 definition_files, config):
        """Construct a namelist generator.

        Arguments:
        `case`             - `Case` object corresponding to the current case.
        `infiles`          - List of files with user namelist options.
        `definition_files` - List of XML files containing namelist definitions.
        `config`           - A dictionary of attributes for matching defaults.
        """
        # Save off important information from inputs.
        self._case = case
        self._din_loc_root = case.get_value('DIN_LOC_ROOT')

        # Create definition object.
        self._definition = NamelistDefinition(definition_files[0], attributes=config)

        # Determine array of _stream_variables from definition object
        # This is only applicable to data models
        self._streams_namelists = {"streams": []}
        self._streams_variables = self._definition.get_per_stream_entries()
        for variable in self._streams_variables:
            self._streams_namelists[variable] = []

        # Create namelist object.
        self._namelist = Namelist()
        for file_ in infiles:
            # Parse settings in "groupless" mode.
            nml_dict = parse(in_file=file_, groupless=True)

            # Add groups using the namelist definition.
            new_namelist = self._definition.dict_to_namelist(nml_dict, filename=file_)

            # Make sure that the input is actually valid.
            self._definition.validate(new_namelist, filename=file_)

            # Merge into existing settings (earlier settings have precedence
            # over later settings).
            self._namelist.merge_nl(new_namelist)

    # Define __enter__ and __exit__ so that we can use this as a context manager
    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def get_definition_nodes(self, skip_groups=[]):
        return self._definition.get_entry_nodes(skip_groups=skip_groups)

    def get_definition_entries(self, skip_groups=[]):
        """Return array of names of all definition entries"""
        return self._definition.get_entries(skip_groups=skip_groups)

    @staticmethod
    def quote_string(string):
        """Convert a string to a quoted Fortran literal.

        Does nothing if the string appears to be quoted already.
        """
        if string == "" or \
           (string[0] not in ('"', "'") or string[0] != string[-1]):
            string = string_to_character_literal(string)
        return string

    def _to_python_value(self, name, literals, node=None):
        """Transform a literal list as needed for `get_value`."""
        var_type, _, var_size, = self._definition.split_type_string(name, self._definition.get_type_info(name, node=node))
        if len(literals) > 0:
            value = expand_literal_list(literals)
        else:
            value = ''
            return value
        for i, scalar in enumerate(value):
            if scalar == '':
                value[i] = None
            elif var_type == 'character':
                value[i] = character_literal_to_string(scalar)
        if var_size == 1:
            return value[0]
        else:
            return value

    def _to_namelist_literals(self, name, value, node=None):
        """Transform a literal list as needed for `set_value`.

        This is the inverse of `_to_python_value`, except that many of the
        changes have potentially already been performed.
        """
        var_type, _, var_size, =  self._definition.split_type_string(name, self._definition.get_type_info(name, node=node))
        if var_size == 1 and not isinstance(value, list):
            value = [value]
        for i, scalar in enumerate(value):
            if scalar is None:
                value[i] = ""
            elif var_type == 'character':
                expect(not isinstance(scalar, list), name)
                value[i] = self.quote_string(scalar)
        return compress_literal_list(value)

    def get_value(self, name, node=None):
        """Get the current value of a given namelist variable.

        Note that the return value of this function is always a string or a list
        of strings. E.g. the scalar logical value .false. will be returned as
        `".false."`, while an array of two .false. values will be returned as
        `[".false.", ".false."]`. Whether or not a value is scalar is determined
        by checking the array size in the namelist definition file.

        Null values are converted to `None`, and repeated values are expanded,
        e.g. `['2*3']` is converted to `['3', '3', '3']`.

        For character variables, the value is converted to a Python string (e.g.
        quotation marks are removed).

        All other literals are returned as the raw string values that will be
        written to the namelist.
        """
        return self._to_python_value(name, self._namelist.get_value(name), node=node)

    def set_value(self, name, value, node=None):
        """Set the current value of a given namelist variable.

        Usually, you should use `add_default` instead of this function.

        The `name` argument is the name of the variable to set, and the `value`
        is a list of strings to use as settings. If the variable is scalar, the
        list is optional; i.e. a scalar logical can be set using either
        `value='.false.'` or `value=['.false.']`. If the variable is of type
        character, and the input is missing quotes, quotes will be added
        automatically. If `None` is provided in place of a string, this will be
        translated to a null value.

        Note that this function will overwrite the current value, which may hold
        a user-specified setting. Even if `value` is (or contains) a null value,
        the old setting for the variable will be thrown out completely.
        """
        if node is not None:
            var_group = self._definition._get_node_element_info(node, "group")
        else:
            var_group = self._definition.get_node_element_info(name, "group")
        literals = self._to_namelist_literals(name, value, node=node)
        self._namelist.set_variable_value(var_group, name, literals)

    def get_default(self, name, config=None, allow_none=False, node=None):
        """Get the value of a variable from the namelist definition file.

        The `config` argument is passed through to the underlying
        `NamelistDefaults.get_value` call as the `attribute` argument.

        The return value of this function is a list of values that were found in
        the defaults file. If there is no matching default, this function
        returns `None` if `allow_none=True` is passed, otherwise an error is
        raised.

        Note that we perform some translation of the values, since there are a
        few differences between Fortran namelist literals and values in the
        defaults file:
        1) In the defaults file, whitespace is ignored except within strings, so
           the output of this function strips out most whitespace. (This implies
           that commas are the only way to separate array elements in the
           defaults file.)
        2) In the defaults file, quotes around character literals (strings) are
           optional, as long as the literal does not contain whitespace, commas,
           or (single or double) quotes. If a setting for a character variable
           does not seem to have quotes (and is not a null value), this function
           will add them.
        3) Default values may refer to variables in a case's `env_*.xml` files.
           This function replaces references of the form `$VAR` or `${VAR}` with
           the value of the variable `VAR` in an env file, if that variable
           exists. This behavior is suppressed within single-quoted strings
           (similar to parameter expansion in shell scripts).
        """
        if node is not None:
            # this uses the entry_id.py _get_value_match function
            default = self._definition._get_value_match(node, attributes=config, exact_match=False)
            if default is None: 
                default = ''
            else:
                default =  self._definition._split_defaults_text(default)
        else:
            # this uses the namelist_definition.py get_value_match function
            default = self._definition.get_value_match(name, attributes=config, exact_match=False)

        if default is None:
            expect(allow_none, "No default value found for %s." % name)
            return None
        default = expand_literal_list(default)

        if node is not None:
            var_type,_,_ = self._definition.split_type_string(name, self._definition.get_type_info(name, node=node))
        else:
            var_type,_,_ = self._definition.split_type_string(name, self._definition.get_type_info(name))

        for i, scalar in enumerate(default):
            # Skip single-quoted strings.
            if var_type == 'character' and scalar != '' and \
               scalar[0] == scalar[-1] == "'":
                continue
            match = _var_ref_re.search(scalar)
            while match:
                env_val = self._case.get_value(match.group('name'))
                expect(env_val is not None,
                       "Namelist default for variable %s refers to unknown XML "
                       "variable %s." % (name, match.group('name')))
                scalar = scalar.replace(match.group(0), str(env_val), 1)
                match = _var_ref_re.search(scalar)
            default[i] = scalar

        # Deal with missing quotes.

        if var_type == 'character':
            for i, scalar in enumerate(default):
                # Preserve null values.
                if scalar != '':
                    default[i] = self.quote_string(scalar)

        default = self._to_python_value(name, default, node=node)
        return default

    def get_streams(self):
        """Get a list of all streams used for the current data model  mode."""
        return self.get_default("streamslist")

    def _sub_fields(self, varnames):
        """Substitute indicators with given values in a list of fields.

        Replace any instance of the following substring indicators with the
        appropriate values:
            %glc = two-digit GLC elevation class from 01 through glc_nec

        The difference between this function and `_sub_paths` is that this
        function is intended to be used for variable names (especially from the
        `strm_datvar` defaults), whereas `_sub_paths` is intended for use on
        input data file paths.

        Returns a string.

        Example: If `_sub_fields` is called with an array containing two
        elements, each of which contains two strings, and glc_nec=3:
             foo               bar
             s2x_Ss_tsrf%glc   tsrf%glc
         then the returned array will be:
             foo               bar
             s2x_Ss_tsrf01     tsrf01
             s2x_Ss_tsrf02     tsrf02
             s2x_Ss_tsrf03     tsrf03
        """
        lines = varnames.split("\n")
        new_lines = []
        for line in lines:
            if not line:
                continue
            if "%glc" in line:
                glc_nec_indices = range(self._case.get_value('GLC_NEC'))
                glc_nec_indices.append(glc_nec_indices[-1] + 1)
                glc_nec_indices.pop(0)
                for i in glc_nec_indices:
                    new_lines.append(line.replace("%glc", "%02d" % i))
            else:
                new_lines.append(line)
        return "\n".join(new_lines)

    @staticmethod
    def _days_in_month(month, year=1):
        """Number of days in the given month (specified as an int, 1-12).

        The `year` argument gives the year for which to request the number of
        days, in a Gregorian calendar. Defaults to `1` (not a leap year).
        """
        month_start = datetime.date(year, month, 1)
        if month == 12:
            next_year = year+1
            next_month = 1
        else:
            next_year = year
            next_month = month + 1
        next_month_start = datetime.date(next_year, next_month, 1)
        return (next_month_start - month_start).days

    def _sub_paths(self, filenames, year_start, year_end):
        """Substitute indicators with given values in a list of filenames.

        Replace any instance of the following substring indicators with the
        appropriate values:
            %y    = year from the range year_start to year_end
            %ym   = year-month from the range year_start to year_end with all 12
                    months
            %ymd  = year-month-day from the range year_start to year_end with
                    all 12 months

        For the date indicators, the year may be prefixed with a number of
        digits to use (the default is 4). E.g. `%2ymd` can be used to change the
        number of year digits from 4 to 2.

        Note that we assume that there is no mixing and matching of date
        indicators, i.e. you cannot use `%4ymd` and `%2y` in the same line. Note
        also that we use a no-leap calendar, i.e. every month has the same
        number of days every year.

        The difference between this function and `_sub_fields` is that this
        function is intended to be used for file names (especially from the
        `strm_datfil` defaults), whereas `_sub_fields` is intended for use on
        variable names.

        Returns a string (filenames separated by newlines).
        """
        lines = [line for line in filenames.split("\n") if line]
        new_lines = []
        for line in lines:
            match = _ymd_re.search(filenames)
            if match is None:
                new_lines.append(line)
                continue
            if match.group('digits'):
                year_format = "%0"+match.group('digits')+"d"
            else:
                year_format = "%04d"
            for year in range(year_start, year_end+1):
                if match.group('day'):
                    for month in range(1, 13):
                        days = self._days_in_month(month)
                        for day in range(1, days+1):
                            date_string = (year_format + "-%02d-%02d") % \
                                          (year, month, day)
                            new_line = line.replace(match.group(0), date_string)
                            new_lines.append(new_line)
                elif match.group('month'):
                    for month in range(1, 13):
                        date_string = (year_format + "-%02d") % (year, month)
                        new_line = line.replace(match.group(0), date_string)
                        new_lines.append(new_line)
                else:
                    date_string = year_format % year
                    new_line = line.replace(match.group(0), date_string)
                    new_lines.append(new_line)
        return "\n".join(new_lines)

    def create_stream_file_and_update_shr_strdata_nml(self, config, #pylint:disable=too-many-locals
                           stream, stream_path, data_list_path):
        """Write the pseudo-XML file corresponding to a given stream.

        Arguments:
        `config` - Used to look up namelist defaults. This is used *in addition*
                   to the `config` used to construct the namelist generator. The
                   main reason to supply additional configuration options here
                   is to specify stream-specific settings.
        `stream` - Name of the stream.
        `stream_path` - Path to write the stream file to.
        `data_list_path` - Path of file to append input data information to.
        """
        # Stream-specific configuration.
        config = config.copy()
        config["stream"] = stream

        # Figure out the details of this stream.
        if stream in ("prescribed", "copyall"):
            # Assume only one file for prescribed mode!
            grid_file = self.get_default("strm_grid_file", config)
            domain_filepath, domain_filenames = os.path.split(grid_file)
            data_file = self.get_default("strm_data_file", config)
            data_filepath, data_filenames = os.path.split(data_file)
        else:
            domain_filepath = self.get_default("strm_domdir", config)
            domain_filenames = self.get_default("strm_domfil", config)
            data_filepath = self.get_default("strm_datdir", config)
            data_filenames = self.get_default("strm_datfil", config)

        domain_varnames = self._sub_fields(self.get_default("strm_domvar", config))
        data_varnames = self._sub_fields(self.get_default("strm_datvar", config))
        offset = self.get_default("strm_offset", config)
        year_start = int(self.get_default("strm_year_start", config))
        year_end = int(self.get_default("strm_year_end", config))
        data_filenames = self._sub_paths(data_filenames, year_start, year_end)
        domain_filenames = self._sub_paths(domain_filenames, year_start, year_end)

        # Overwrite domain_file if should be set from stream data
        if domain_filenames == 'null':
            domain_filepath = data_filepath
            domain_filenames = data_filenames.splitlines()[0]

        stream_file_text = _stream_file_template.format(
            domain_varnames=domain_varnames,
            domain_filepath=domain_filepath,
            domain_filenames=domain_filenames,
            data_varnames=data_varnames,
            data_filepath=data_filepath,
            data_filenames=data_filenames,
            offset=offset,
        )

        with open(stream_path, 'w') as stream_file:
            stream_file.write(stream_file_text)
        with open(data_list_path, 'a') as input_data_list:
            for i, filename in enumerate(domain_filenames.split("\n")):
                if filename.strip() == '':
                    continue
                filepath = os.path.join(domain_filepath, filename.strip())
                input_data_list.write("domain%d = %s\n" % (i+1, filepath))
            for i, filename in enumerate(data_filenames.split("\n")):
                if filename.strip() == '':
                    continue
                filepath = os.path.join(data_filepath, filename.strip())
                input_data_list.write("file%d = %s\n" % (i+1, filepath))
        self.update_shr_strdata_nml(config, stream, stream_path)

    def update_shr_strdata_nml(self, config, stream, stream_path):
        """Updates values for the `shr_strdata_nml` namelist group.

        This should be done once per stream, and it shouldn't usually be called
        directly, since `create_stream_file` calls this method itself.
        """
        assert config['stream'] == stream, \
            "config stream is %s, but input stream is %s" % \
            (config['stream'], stream)
        # Double-check the years for sanity.
        year_start = int(self.get_default("strm_year_start", config))
        year_end = int(self.get_default("strm_year_end", config))
        year_align = int(self.get_default("strm_year_align", config))
        expect(year_end >= year_start,
               "Stream %s starts at year %d, but ends at earlier year %d." %
               (stream, year_start, year_end))
        # Add to streams file.
        stream_string = "%s %d %d %d" % (os.path.basename(stream_path),
                                         year_align, year_start, year_end)
        self._streams_namelists["streams"].append(stream_string)
        for variable in self._streams_variables:
            default = self.get_default(variable, config)
            expect(len(default) == 1,
                   "Stream %s had multiple settings for variable %s." %
                   (stream, variable))
            self._streams_namelists[variable].append(default[0])

    def set_abs_file_path(self, file_path):
        """If `file_path` is relative, make it absolute using `DIN_LOC_ROOT`.

        If an absolute path is input, it is returned unchanged.
        """
        if os.path.isabs(file_path):
            return file_path
        else:
            fullpath = os.path.join(self._din_loc_root, file_path)
            return fullpath

    def add_default(self, name, value=None, node=None):
        """Add a value for the specified variable to the namelist.

        If the specified variable is already defined in the object, the existing
        value is preserved. Otherwise, the `value` argument, if provided, will
        be used to set the value. If no such value is found, the defaults file
        will be consulted. If null values are present in any of the above, the
        result will be a merged array of values.

        If no value for the variable is found via any of the above, this method
        will raise an exception.
        """
        if node is None:
            nodes = self._definition.get_nodes_by_id(name)
            expect(len(nodes) == 1, "incorrect number of matchs %s"%len(nodes))
            node = nodes[0]

        group = self._definition._get_node_element_info(node, "group")

        # Use this to see if we need to raise an error when nothing is found.
        have_value = False

        # Check for existing value.
        current_literals = self._namelist.get_variable_value(group, name)
        if current_literals != [""]:
            have_value = True
            # Do not proceed further since this has been obtained the -infile contents
            return 

        # Check for input argument.
        if value is not None:
            have_value = True
            # if compression were to occur, this is where it does
            literals = self._to_namelist_literals(name, value) 
            current_literals = merge_literal_lists(literals, current_literals)

        # Check for default value.
        default = self.get_default(name, allow_none=True, node=node)
        if default is not None:
            have_value = True
            default_literals = self._to_namelist_literals(name, default)
            current_literals = merge_literal_lists(default_literals, current_literals)
        expect(have_value, "No default value found for %s." % name)

        # Go through file names and prepend input data root directory for
        # absolute pathnames.
        var_input_pathname = self._definition.get_input_pathname(name)
        if var_input_pathname == 'abs':
            current_literals = expand_literal_list(current_literals)
            for i, literal in enumerate(current_literals):
                if literal == '':
                    continue
                file_path = character_literal_to_string(literal)
                # NOTE - these are hard-coded here and a better way is to make these extensible
                if file_path == 'UNSET' or file_path == 'idmap':
                    continue
                if file_path == 'null':
                    continue
                file_path = self.set_abs_file_path(file_path)
                if not os.path.exists(file_path):
                    logger.warn ("File not found: %s = %s, will attempt to download in check_input_data phase" % (name, literal))
                current_literals[i] = string_to_character_literal(file_path)
            current_literals = compress_literal_list(current_literals)

        # Set the new value.
        self._namelist.set_variable_value(group, name, current_literals)

    def create_shr_strdata_nml(self):
        """Set defaults for `shr_strdata_nml` variables other than the variable domainfile """
        self.add_default("datamode")
        if self.get_value("datamode") != 'NULL':
            self.add_default("streams",
                             value=self._streams_namelists["streams"])
            for variable in self._streams_variables:
                self.add_default(variable,
                                 value=self._streams_namelists[variable])

    def _write_input_files(self, input_data_list):
        """Write input data files to list."""
        for group_name in self._namelist.get_group_names():
            for variable_name in self._namelist.get_variable_names(group_name):
                input_pathname = self._definition.get_node_element_info(variable_name, "input_pathname")
                if input_pathname is not None:
                    # This is where we end up for all variables that are paths
                    # to input data files.
                    literals = self._namelist.get_variable_value(group_name,
                                                                 variable_name)
                    for literal in literals:
                        file_path = character_literal_to_string(literal)
                        # NOTE - these are hard-coded here and a better way is to make these extensible
                        if file_path == 'UNSET' or file_path == 'idmap':
                            continue
                        if input_pathname == 'abs':
                            # No further mangling needed for absolute paths.
                            pass
                        elif input_pathname.startswith('rel:'):
                            # The part past "rel" is the name of a variable that
                            # this variable specifies its path relative to.
                            root_var = input_pathname[4:]
                            root_dir = self.get_value(root_var)
                            file_path = os.path.join(root_dir, file_path)
                        else:
                            expect(False,
                                   "Bad input_pathname value: %s." %
                                   input_pathname)
                        # Write to the input data list.
                        input_data_list.write("%s = %s\n" % (variable_name, file_path))

    def write_output_file(self, namelist_file, data_list_path=None, groups=None, sorted_groups=True):
        """Write out the namelists and input data files.

        The `namelist_file` and `modelio_file` are the locations to which the
        component and modelio namelists will be written, respectively. The
        `data_list_path` argument is the location of the `*.input_data_list`
        file, which will have the input data files added to it.
        """
        self._definition.validate(self._namelist)
        if groups is None:
            groups = self._namelist.get_group_names()

        # remove groups that are never in namelist file
        if "modelio" in groups:
            groups.remove("modelio")
        if "seq_maps" in groups:
            groups.remove("seq_maps")

        # write namelist file
        self._namelist.write(namelist_file, groups=groups, sorted_groups=sorted_groups)

        if data_list_path is not None:
            # append to input_data_list file
            with open(data_list_path, "a") as input_data_list:
                self._write_input_files(input_data_list)

    def add_nmlcontents(self, filename, group, append=True, format_="nmlcontents", sorted_groups=True):
        """ Write only contents of nml group """
        self._namelist.write(filename, groups=[group], append=append, format_=format_, sorted_groups=sorted_groups)

    def write_seq_maps(self, filename):
        """ Write out seq_maps.rc"""
        self._namelist.write(filename, groups=["seq_maps"], format_="rc")

    def write_modelio_file(self, filename):
        """ Write  component modelio files"""
        self._namelist.write(filename, groups=["modelio", "pio_inparm"], format_="nml")
