"""
Classes used to build the CIME Macros file.

The main "public" class here is Build. It is initialized with machine-specific
information, and its write_macros method is the driver for translating the
config_build.xml file into a Makefile or CMake-format Macros file.

For developers, here's the role of the other classes in the process:

- A CompilerBlock is responsible for translating the XML code in a <compiler>
  tag into Python data structures.

- A PossibleValues object keeps track of all the settings that could affect a
  particular variable, and is the main way that these settings are stored.

- A MacroConditionTree is the structure that is responsible for writing out the
  settings. While the PossibleValues objects are organized by variable name, the
  MacroConditionTree is organized by conditional blocks, and thus roughly
  plays the role of a syntax tree corresponding to the Makefile/CMake output.

In more detail:

- Build.write_macros immediately creates a MakeMacroWriter or CMakeMacroWriter
  to translate strings for the build system.

- It also creates value_lists, a dictionary of PossibleValues objects, with
  variable names as the keys. Each variable has a single PossibleValues object
  associated with it.

- For each <compiler> element, Build.write_macros creates a CompilerBlock
  instance. This object is responsible for translating the XML in its block, in
  order to populate the PossibleValues instances. This includes handling the
  <var>/<env>/<shell> tags, and keeping track of dependencies induced by one
  variable referencing another's value.

- The PossibleValues object holds the information about how one variable can be
  set, based on various build options. It has two main roles:
   1. As we iterate through the XML input file, each setting is added to the
      relevant PossibleValues object. The PossibleValues object contains lists
      of settings sorted by how machine-specific those settings are.
   2. The PossibleValues object iterates through the list of settings to check
      for ambiguities. E.g. if there is a setting for DEBUG=TRUE, and another
      setting for MPILIB=mpi-serial, it is ambiguous in the case where both
      conditions hold.

- A ValueSetting object is a simple struct that a setting from the XML file is
  translated to. The lists in the PossibleValues class contain these objects.

- Once the XML has all been read in and the PossibleValues objects are
  populated, the dependencies among variables are checked in Build.write_macros.
  For each variable, if all its dependencies have been handled, it is converted
  to a MacroConditionTree merged with all other trees for variables that are
  ready, and written out. Then we loop through the variable list again to check
  for variables whose dependencies are all handled.

- The MacroConditionTree acts as a primitive syntax tree. Its __init__ method
  reorganizes the data into conditional blocks, and its write_out method writes
  uses the MakeMacroWriter/CMakeMacroWrite object to write to the Macros file.
  MacroConditionTree objects can be merged to reduce the length of the output.
"""

# These don't seem to be particularly useful checks.
# pylint: disable=invalid-name,too-few-public-methods,unused-wildcard-import
# pylint: disable=wildcard-import

from CIME.macros_writers import *
from CIME.utils import get_cime_root
from CIME.XML.machines import Machines # pylint: disable=unused-import
from CIME.XML.standard_module_setup import *

__all__ = ["Build"]

logger = logging.getLogger(__name__)

class ValueSetting(object):

    """Holds data about how a value can be assigned to a variable.

    Note that this class doesn't know or care *which* variable might be
    assigned in this way, only that there is a procedure to perform that
    operation

    Public attributes:
    value - The actual value that will be set.
    do_append - Boolean describing whether the value should be
                appended to the existing value of the variable rather
                than overwriting other settings.
    conditions - Dictionary containing the set of values that different
                 variables have to have to use this setting (e.g.
                 DEBUG="TRUE" might be a condition on a debug flag).
    set_up - List of any commands that have to be executed in the build
             system before this setting can occur.
    tear_down - List of any commands that should be executed to clean up
                after setting the variable.

    Public methods:
    is_ambiguous_with
    """

    def __init__(self, value, do_append, conditions, set_up, tear_down): #  pylint: disable=too-many-arguments
        """Create a ValueSetting object by specifying all its data."""
        self.value = value
        self.do_append = do_append
        self.conditions = conditions
        self.set_up = set_up
        self.tear_down = tear_down

    def is_ambiguous_with(self, other):
        """Check to see if this setting conflicts with another one.

        The purpose of this routine is to see if two settings can coexist
        in the same Macros file, or if doing so would raise an ambiguity
        about which one should be preferred over the other. Note that this
        is a symmetric relation (this function returns the same value if
        self and other are swapped).

        The rules to determine this are as follows:

        1) If one or both settings are appending to the value, there's no
           ambiguity, because both can cooperate to set the value.

        >>> a = ValueSetting('foo', True, dict(), [], [])
        >>> b = ValueSetting('bar', False, dict(), [], [])
        >>> a.is_ambiguous_with(b)
        False
        >>> b.is_ambiguous_with(a)
        False

        2) If the two settings have conflicting conditions, then there
           is no ambiguity because they can't both apply to the same
           build.

        >>> a = ValueSetting('foo', False, {"DEBUG": "TRUE"}, [], [])
        >>> b = ValueSetting('bar', False, {"DEBUG": "FALSE"}, [], [])
        >>> a.is_ambiguous_with(b)
        False

        3) If one setting is strictly more specific than the other, then
           there's no ambiguity, because we prefer the more specific
           setting whenever both apply to a build.

        >>> a = ValueSetting('foo', False, {"DEBUG": "TRUE"}, [], [])
        >>> b = ValueSetting('bar', False, {"DEBUG": "TRUE", "MPILIB": "mpich2"}, [], [])
        >>> a.is_ambiguous_with(b)
        False
        >>> b.is_ambiguous_with(a)
        False

        4) All other cases are considered ambiguous.

        >>> a = ValueSetting('foo', False, dict(), [], [])
        >>> b = ValueSetting('bar', False, dict(), [], [])
        >>> a.is_ambiguous_with(b)
        True
        >>> a = ValueSetting('foo', False, {"DEBUG": "TRUE"}, [], [])
        >>> b = ValueSetting('bar', False, {"MPILIB": "mpich2"}, [], [])
        >>> a.is_ambiguous_with(b)
        True
        """
        # Append check.
        if self.do_append or other.do_append:
            return False
        # Consistency check.
        for var_name in self.conditions:
            if var_name not in other.conditions:
                continue
            if self.conditions[var_name] != other.conditions[var_name]:
                return False
        # Specificity check.
        # One setting being more specific than the other is equivalent to
        # its set of conditions being a proper superset of the others.
        self_set = set(self.conditions.keys())
        other_set = set(other.conditions.keys())
        if self_set < other_set or self_set > other_set:
            return False
        # Any situation we couldn't resolve is ambiguous.
        return True


class PossibleValues(object):

    """Holds a list of settings for a single "Macros" variable.

    This helper class takes in variable settings and, for each one, decides
    whether to throw it out, add it to the list of values, or replace the
    existing list of values with the new, more specific setting.

    This class also performs ambiguity checking; if it is possible at build
    time for more than one setting to match the same variable, this is
    considered an error.

    Public attributes:
    name - The name of the variable.
    settings - The current list of possible initial settings for the
               variable.
    append_settings - A dictionary of lists of possible appending settings
                      for the variable, with the specificity of each list
                      as the associated dictionary key.
    depends - The current list of variables that this variable depends on
              to get its value.

    Public methods:
    add_setting
    ambiguity_check
    to_cond_trees
    """

    def __init__(self, name, setting, specificity, depends):
        """Construct a PossibleValues object.

        The name argument refers to the name of the variable. The other
        arguments are the same as for append_match.
        """
        self.name = name
        self.depends = depends
        # If this is an appending setting, its specificity can't cause it
        # to overwrite other settings, but we want to keep track of it.
        if setting.do_append:
            self.settings = []
            self.append_settings = {specificity: [setting]}
            self._specificity = 0
        else:
            self.settings = [setting]
            self.append_settings = {}
            self._specificity = specificity

    def add_setting(self, setting, specificity, depends):
        """Add a possible value for a variable.

        Arguments:
        setting - A ValueSetting to start the list.
        specificity - An integer representing how specific the setting is.
                      Only the initial settings with the highest
                      specificity and appending settings with at least that
                      specificity will actually be kept in the list. The
                      lowest allowed specificity is 0.
        depends - A set of variable names, specifying the variables that
                  have to be set before this setting can be used (e.g. if
                  SLIBS refers to NETCDF_PATH, then NETCDF_PATH has to be
                  set first).

        >>> a = ValueSetting('foo', False, dict(), [], [])
        >>> b = ValueSetting('bar', False, dict(), [], [])
        >>> vals = PossibleValues('var', a, 0, {'dep1'})
        >>> vals.add_setting(b, 1, {'dep2'})
        >>> a not in vals.settings and b in vals.settings
        True
        >>> 'dep1' not in vals.depends and 'dep2' in vals.depends
        True
        >>> vals.add_setting(a, 1, {'dep1'})
        >>> a in vals.settings and b in vals.settings
        True
        >>> 'dep1' in vals.depends and 'dep2' in vals.depends
        True
        """
        if setting.do_append:
            # Appending settings with at least the current level of
            # specificity should be kept.
            if specificity >= self._specificity:
                if specificity not in self.append_settings:
                    self.append_settings[specificity] = []
                self.append_settings[specificity].append(setting)
                self.depends |= depends
        else:
            # Add equally specific settings to the list.
            if specificity == self._specificity:
                self.settings.append(setting)
                self.depends |= depends
            # Replace the list if the setting is more specific.
            elif specificity > self._specificity:
                self.settings = [setting]
                self._specificity = specificity
                self.depends = depends
        # Do nothing if the setting is less specific.

    def ambiguity_check(self):
        """Check the current list of settings for ambiguity.

        This function raises an error if an ambiguity is found.
        """
        for i in range(len(self.settings)-1):
            for other in self.settings[i+1:]:
                expect(not self.settings[i].is_ambiguous_with(other),
                       "Variable "+self.name+" is set ambiguously in "
                       "config_build.xml. Check the file for these "
                       "conflicting settings: \n1: {}\n2: {}".format(
                           self.settings[i].conditions, other.conditions))

    def to_cond_trees(self):
        """Convert this object to a pair of MacroConditionTree objects.

        This represents the step where the list of possible values is
        frozen and we're ready to convert it into an actual text file. This
        object is checked for ambiguities before conversion.

        The return value is a tuple of two trees. The first contains all
        initial settings, and the second contains all appending settings.
        If either would be empty, None is returned instead.
        """
        self.ambiguity_check()
        if self.settings:
            normal_tree = MacroConditionTree(self.name, self.settings)
        else:
            normal_tree = None
        append_settings = []
        for specificity in self.append_settings:
            if specificity >= self._specificity:
                append_settings += self.append_settings[specificity]
        if append_settings:
            append_tree = MacroConditionTree(self.name, append_settings)
        else:
            append_tree = None
        return (normal_tree, append_tree)


class MacroConditionTree(object): # pylint: disable=too-many-instance-attributes

    """Tree containing the various possible settings of a specific macro.

    Unlike the PossibleValues class, this class assumes that we have
    finished determining which settings could apply on a given machine. It
    also sorts the settings based on the conditions under which they take
    effect, in preparation for writing out the Macros file itself.

    Public methods:
    merge
    write_out
    """

    def __init__(self, name, settings):
        """Create a MacroConditionTree recursively.

        Arguments:
        name - Name of the variable.
        settings - A list of all settings for this variable.
        """
        # Search for any conditions controlling the number of settings.
        condition = None
        # Prefer the COMPILER attribute as the top level attribute, for
        # readability of the merged file.
        if any("COMPILER" in setting.conditions for setting in settings):
            condition = "COMPILER"
        else:
            # To make merging more effective, sort the conditions.
            all_conditions = []
            for setting in settings:
                all_conditions += setting.conditions.keys()
            if all_conditions:
                condition = sorted(all_conditions)[0]
        if condition is None:
            # If there are no conditions, we have reached a leaf.
            # We combine whatever settings are left; there should be at
            # most one non-appending setting, or an arbitrary number of
            # appending settings.
            self._is_leaf = True
            self._assignments = []
            self._set_up = []
            self._tear_down = []
            self._do_append = True
            for setting in settings:
                if not setting.do_append:
                    self._do_append = False
                    assert len(settings) == 1, \
                        "Internal error in macros: An ambiguity was " \
                        "found after the ambiguity check was complete, " \
                        "or there is a mixture of appending and initial " \
                        "settings in the condition tree."
                self._assignments.append((name, setting.value))
                self._set_up += setting.set_up
                self._tear_down += setting.tear_down
        else:
            # If a condition was found, partition the settings depending on
            # how they use it, and recursively create a tree for each
            # partition.
            self._is_leaf = False
            self._condition = condition
            partition = dict()
            for setting in settings:
                # If some of the settings don't use a condition, we use
                # None to represent that.
                cond_val = setting.conditions.pop(condition, None)
                if cond_val in partition:
                    partition[cond_val].append(setting)
                else:
                    partition[cond_val] = [setting]
            branches = dict()
            for cond_val in partition:
                branches[cond_val] = \
                            MacroConditionTree(name, partition[cond_val])
            self._branches = branches

    # pylint shouldn't concern itself with the way that we access other, since
    # it's actually a member of the same class.
    # pylint:disable=protected-access
    def merge(self, other):
        """Merge another tree with this one.

        This should be considered destructive to both trees. The only valid
        value is the one that's returned.
        """
        if self._is_leaf:
            if other._is_leaf:
                assert self._do_append == other._do_append, \
                    "Internal error in macros: Tried to merge an " \
                    "appending tree with a tree containing initial "\
                    "settings."
                # If both are leaves, just merge the values.
                self._assignments += other._assignments
                self._set_up += other._set_up
                self._tear_down += other._tear_down
                return self
            else:
                # If other is not a leaf, swap the arguments so that self
                # is the one that's not a leaf, handled below.
                return other.merge(self)
        else:
            # If self is not a leaf but other is, it should go in
            # self._branches[None]. The same goes for the case where the
            # conditions don't match, and self._condition is last
            # alphabetically.
            if other._is_leaf or self._condition > other._condition:
                if None in self._branches:
                    self._branches[None] = self._branches[None].merge(other)
                else:
                    self._branches[None] = other
                return self
            else:
                # If the other condition comes last alphabetically, swap
                # the order.
                if self._condition < other._condition:
                    return other.merge(self)
                # If neither is a leaf and their conditions match, merge
                # their sets of branches.
                for (cond_val, other_branch) in other._branches.items():
                    if cond_val in self._branches:
                        self._branches[cond_val] = \
                            self._branches[cond_val].merge(other_branch)
                    else:
                        self._branches[cond_val] = other_branch
                return self
    # pylint:enable=protected-access

    def write_out(self, writer):
        """Write tree to file.

        The writer argument is an object inheriting from MacroWriterBase.
        This function first writes out all the initial settings with
        appropriate conditionals, then the appending settings.
        """
        if self._is_leaf:
            for line in self._set_up:
                writer.write_line(line)
            for (name, value) in self._assignments:
                if self._do_append:
                    writer.append_variable(name, value)
                else:
                    writer.set_variable(name, value)
            for line in self._tear_down:
                writer.write_line(line)
        else:
            condition = self._condition
            # Take care of the settings that don't use this condition.
            if None in self._branches:
                self._branches[None].write_out(writer)
            # Now all the if statements.
            for cond_val in self._branches:
                if cond_val is None:
                    continue
                env_ref = writer.environment_variable_string(condition)
                writer.start_ifeq(env_ref, cond_val)
                self._branches[cond_val].write_out(writer)
                writer.end_ifeq()


def _merge_optional_trees(tree, big_tree):
    """Merge two MacroConditionTrees when one or both objects may be `None`."""
    if tree is not None:
        if big_tree is None:
            return tree
        else:
            return big_tree.merge(tree)
    else:
        return big_tree


class CompilerBlock(object):

    """Data used to translate a single <compiler> element.

    This is used during write_macros to traverse the XML and create a list
    of settings specified in the element.

    Public methods:
    add_settings_to_lists
    matches_machine
    """

    def __init__(self, writer, compiler_elem, machobj):
        """Construct a CompilerBlock.

        Arguments:
        writer - The Makefile/CMake writer object.
        compiler_elem - An xml.ElementTree.Element corresponding to this
                        <compiler> element.
        machobj - Machines object for this machine.
        """
        self._writer = writer
        self._compiler_elem = compiler_elem
        self._machobj = machobj
        # If there's no COMPILER attribute, self._compiler is None.
        self._compiler = compiler_elem.get("COMPILER")
        self._specificity = 0

    def _handle_references(self, elem, set_up, tear_down, depends):
        """Expand markup used internally.

        This function is responsible for expanding <env>, <var>, and
        <shell> tags into Makefile/CMake syntax.

        Arguments:
        elem - An ElementTree.Element containing text to expand.
        set_up - A list to add any preparation commands to.
        tear_down - A list to add any cleanup commands to.
        depends - A set of variables that need to be set before this one.

        Note that while the return value of this function is the expanded
        text, the set_up, tear_down, and depends variables are also
        modified and thus serve as additional outputs.
        """
        writer = self._writer
        output = elem.text
        if output is None:
            output = ""
        for child in elem:
            if child.tag == "env":
                # <env> tags just need to be expanded by the writer.
                output += writer.environment_variable_string(child.text)
            elif child.tag == "shell":
                # <shell> tags can contain other tags, so handle those.
                command = self._handle_references(child, set_up, tear_down,
                                                  depends)
                new_set_up, inline, new_tear_down = \
                                    writer.shell_command_strings(command)
                output += inline
                if new_set_up is not None:
                    set_up.append(new_set_up)
                if new_tear_down is not None:
                    tear_down.append(new_tear_down)
            elif child.tag == "var":
                # <var> commands also need expansion by the writer, and can
                # add dependencies.
                var_name = child.text
                output += writer.variable_string(var_name)
                depends.add(var_name)
            else:
                expect(False,
                       "Unexpected tag "+child.tag+" encountered in "
                       "config_build.xml. Check that the file is valid "
                       "according to the schema.")
            if child.tail is not None:
                output += child.tail
        return output

    def _elem_to_setting(self, elem):
        """Take an element and convert it to a ValueSetting.

        Arguments:
        elem - An ElementTree.Element with data to add.

        This function returns a tuple containing a ValueSetting
        corresponding to the element, along with a set of names of
        variables that this setting depends on.
        """
        # Attributes on an element are the conditions on that element.
        conditions = dict(elem.items())
        if self._compiler is not None:
            conditions["COMPILER"] = self._compiler
        # Deal with internal markup.
        set_up = []
        tear_down = []
        depends = set()
        value_text = self._handle_references(elem, set_up,
                                             tear_down, depends)
        # Create the setting object.
        setting = ValueSetting(value_text, elem.tag == "append",
                               conditions, set_up, tear_down)
        return (setting, depends)

    def _add_elem_to_lists(self, name, elem, value_lists):
        """Add an element's data to an appropriate list of value settings.

        Arguments:
        name - The name of the variable being set by this element.
        elem - The element to translate into a ValueSetting.
        value_lists - A dictionary of PossibleValues, containing the lists
                      of all settings for each variable.
        """
        # Skip this if the element's MPILIB is not valid.
        if "MPILIB" in elem.keys() and \
           not self._machobj.is_valid_MPIlib(elem.get("MPILIB")):
            return
        setting, depends = self._elem_to_setting(elem)
        if name not in value_lists:
            value_lists[name] = PossibleValues(name, setting,
                                               self._specificity, depends)
        else:
            value_lists[name].add_setting(setting, self._specificity,
                                          depends)

    def add_settings_to_lists(self, flag_vars, value_lists):
        """Add all data in the <compiler> element to lists of settings.

        Arguments:
        flag_vars - A set of variables containing "flag-like" data.
        value_lists - A dictionary of PossibleValues, containing the lists
                      of all settings for each variable.
        """
        for elem in self._compiler_elem:
            # Deal with "flag"-type variables.
            if elem.tag in flag_vars:
                for child in elem:
                    self._add_elem_to_lists(elem.tag, child, value_lists)
            else:
                self._add_elem_to_lists(elem.tag, elem, value_lists)

    def matches_machine(self):
        """Check whether this block matches a machine/os.

        This also sets the specificity of the block, so this must be called
        before add_settings_to_lists if machine-specific output is needed.
        """
        self._specificity = 0
        if "MACH" in self._compiler_elem.keys():
            if self._machobj.get_machine_name() == \
               self._compiler_elem.get("MACH"):
                self._specificity += 2
            else:
                return False
        if "OS" in self._compiler_elem.keys():
            if self._machobj.get_value("OS") == self._compiler_elem.get("OS"):
                self._specificity += 1
            else:
                return False
        # Check if the compiler is valid on this machine.
        if self._compiler is not None:
            return self._machobj.is_valid_compiler(self._compiler)
        else:
            return True


class Build(object):

    """Class to convert config_build.xml input into a macros file.

    Public attributes:
    os - Operating system used in config_build lookup and ranking.
    machobj - Machines object used in config_build lookup.
    flag_vars - A set of all variables in config_build that contain "flag-
                like" data (i.e. a space-separated list of arguments).

    Public methods:
    write_macros
    """

    def __init__(self, machobj, schema_path=None):
        """Construct a Build given machine-specific information.

        In the process some information about possible variables is read in
        from the schema file.

        Arguments:
        machobj - A Machines object for this machine.
        schema_path (optional) - Path to config_build.xsd within CIME.

        >>> "CFLAGS" in Build('MyMach').flag_vars
        True
        >>> "MPICC" in Build('MyMach').flag_vars
        False
        """
        self.machobj = machobj

        # The schema is used to figure out which variables contain
        # command-line arguments (e.g. compiler flags), since these are
        # processed in a more complex manner than other variables.
        if schema_path is None:
            schema_path = os.path.join(get_cime_root(), "cime_config",
                                       "xml_schemas", "config_build.xsd")

        # Run an XPath query to extract the list of flag variable names.
        ns = {"xs": "http://www.w3.org/2001/XMLSchema"}
        flag_xpath = ".//xs:group[@name='compilerVars']/xs:choice/xs:element[@type='flagsVar']"
        flag_elems = ET.parse(schema_path).getroot().findall(flag_xpath, ns)
        self.flag_vars = set(elem.get('name') for elem in flag_elems)

    def write_macros(self, build_system, xml_file, output):
        """Write a Macros file for this machine.

        Arguments:
        build_system - Format of the file to be written. Currently the only
                       valid values are "Makefile" and "CMake".
        xml_file - File name or file object containing the config_build
                   specification of machine-specific settings.
        output - Text I/O object (inheriting from io.TextIOBase) that
                 output should be written to. Typically, this will be the
                 Macros file, opened for writing.
        """
        # Set up writer for this build system.
        if build_system == "Makefile":
            writer = MakeMacroWriter(output)
        elif build_system == "CMake":
            writer = CMakeMacroWriter(output)
        else:
            expect(False,
                   "Unrecognized build system provided to write_macros: " +
                   build_system)

        # Start processing the file.
        value_lists = dict()
        for compiler_elem in ET.parse(xml_file).findall("compiler"):
            block = CompilerBlock(writer, compiler_elem, self.machobj)
            # If this block matches machine settings, use it.
            if block.matches_machine():
                block.add_settings_to_lists(self.flag_vars, value_lists)

        # Now that we've scanned through the input, output the variable
        # settings.
        vars_written = set()
        while value_lists:
            # Variables that are ready to be written.
            ready_variables = [
                var_name for var_name in value_lists.keys()
                if value_lists[var_name].depends <= vars_written
            ]
            expect(len(ready_variables) > 0,
                   "The config_build XML has bad <var> references. "
                   "Check for circular references or variables that "
                   "are in a <var> tag but not actually defined.")
            big_normal_tree = None
            big_append_tree = None
            for var_name in ready_variables:
                # Note that we're writing this variable.
                vars_written.add(var_name)
                # Make the conditional trees and write them out.
                normal_tree, append_tree = \
                    value_lists[var_name].to_cond_trees()
                big_normal_tree = _merge_optional_trees(normal_tree,
                                                        big_normal_tree)
                big_append_tree = _merge_optional_trees(append_tree,
                                                        big_append_tree)
                # Remove this variable from the list of variables to handle
                # next iteration.
                del value_lists[var_name]
            if big_normal_tree is not None:
                big_normal_tree.write_out(writer)
            if big_append_tree is not None:
                big_append_tree.write_out(writer)
