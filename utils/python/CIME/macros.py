"""
Classes used to build the CIME Macros file.

The main "public" class here is MacroMaker. Other classes are mainly for
internal use, though the various "MacroWriter" classes might be useful
elsewhere one day.
"""

# These don't seem to be particularly useful checks.
# pylint: disable=invalid-name,too-few-public-methods,unused-wildcard-import

from abc import ABCMeta, abstractmethod

from CIME.XML.standard_module_setup import * # pylint: disable=wildcard-import
from CIME.utils import get_cime_root

__all__ = ["MacroMaker"]

logger = logging.getLogger(__name__)

class MacroWriterBase(object):

    """Abstract base class for macro file writers.

    Public attributes:
    indent_increment - Number of spaces to indent if blocks (does not apply
                       to format-specific indentation, e.g. cases where
                       Makefiles must use tabs).
    output - File-like object that output is written to.

    Public methods:
    indent
    indent_left
    indent_right
    write_line
    environment_variable_string
    shell_command_string
    variable_string
    set_variable
    append_variable
    start_ifeq
    end_ifeq
    """

    __metaclass__ = ABCMeta

    indent_increment = 2

    def __init__(self, output):
        """Initialize a macro writer.

        Arguments:
        output - File-like object (probably an io.TextIOWrapper), which
                 will be written to.
        """
        self.output = output
        self._indent_num = 0

    def indent(self):
        """Return an appropriate number of spaces for the indent."""
        return ' ' * self._indent_num

    def indent_left(self):
        """Decrease the amount of line indent."""
        self._indent_num -= 2

    def indent_right(self):
        """Increase the amount of line indent."""
        self._indent_num += 2

    def write_line(self, line):
        """Write a single line of output, appropriately indented.

        A trailing newline is added, whether or not the input has one.
        """
        self.output.write(unicode(self.indent() + line + "\n"))

    @abstractmethod
    def environment_variable_string(self, name):
        """Return an environment variable reference."""
        pass

    @abstractmethod
    def shell_command_strings(self, command):
        """Return strings used to get the output of a shell command.

        Implementations should return a tuple of three strings:
        1. A line that is needed to get the output of the command (or None,
           if a command can be run inline).
        2. A string that can be used within a line to refer to the output.
        3. A line that does any cleanup of temporary variables (or None, if
           no cleanup is necessary).

        Example usage:

        # Get strings and write initial command.
        (pre, var, post) = writer.shell_command_strings(command)
        if pre is not None:
            writer.write(pre)

        # Use the variable to write an if block.
        writer.start_ifeq(var, "TRUE")
        writer.set_variable("foo", "bar")
        writer.end_ifeq()

        # Cleanup
        if post is not None:
            writer.write(post)
        """
        pass

    @abstractmethod
    def variable_string(self, name):
        """Return a string to refer to a variable with the given name."""
        pass

    @abstractmethod
    def set_variable(self, name, value):
        """Write out a statement setting a variable to some value."""
        pass

    def append_variable(self, name, value):
        """Write out a statement appending a value to a string variable."""
        var_string = self.variable_string(name)
        self.set_variable(name, var_string + " " + value)

    @abstractmethod
    def start_ifeq(self, left, right):
        """Write out a statement to start a conditional block.

        The arguments to this method are compared, and the block is entered
        only if they are equal.
        """
        pass

    @abstractmethod
    def end_ifeq(self):
        """Write out a statement to end a block started with start_ifeq."""
        pass


class MakeMacroWriter(MacroWriterBase):

    """Macro writer for the Makefile format.

    For details on the provided methods, see MacroWriterBase, which this
    class inherits from.
    """

    def environment_variable_string(self, name):
        """Return an environment variable reference.

        >>> import io
        >>> s = io.StringIO()
        >>> MakeMacroWriter(s).environment_variable_string("foo")
        '$(foo)'
        """
        return "$(" + name + ")"

    def shell_command_strings(self, command):
        """Return strings used to get the output of a shell command.

        >>> import io
        >>> s = io.StringIO()
        >>> MakeMacroWriter(s).shell_command_strings("echo bar")
        (None, '$(shell echo bar)', None)
        """
        return (None, "$(shell " + command + ")", None)

    def variable_string(self, name):
        """Return a string to refer to a variable with the given name.

        >>> import io
        >>> s = io.StringIO()
        >>> MakeMacroWriter(s).variable_string("foo")
        '$(foo)'
        """
        return "$(" + name + ")"

    def set_variable(self, name, value):
        """Write out a statement setting a variable to some value.

        >>> import io
        >>> s = io.StringIO()
        >>> MakeMacroWriter(s).set_variable("foo", "bar")
        >>> s.getvalue()
        u'foo := bar\\n'
        """
        # Note that ":=" is used so that we can control the behavior for
        # both Makefile and CMake variables similarly.
        self.write_line(name + " := " + value)

    def start_ifeq(self, left, right):
        """Write out a statement to start a conditional block.

        >>> import io
        >>> s = io.StringIO()
        >>> MakeMacroWriter(s).start_ifeq("foo", "bar")
        >>> s.getvalue()
        u'ifeq (foo,bar)\\n'
        """
        self.write_line("ifeq (" + left + "," + right + ")")
        self.indent_right()

    def end_ifeq(self):
        """Write out a statement to end a block started with start_ifeq.

        >>> import io
        >>> s = io.StringIO()
        >>> writer = MakeMacroWriter(s)
        >>> writer.start_ifeq("foo", "bar")
        >>> writer.set_variable("foo2", "bar2")
        >>> writer.end_ifeq()
        >>> s.getvalue()
        u'ifeq (foo,bar)\\n  foo2 := bar2\\nendif\\n'
        """
        self.indent_left()
        self.write_line("endif")


class CMakeMacroWriter(MacroWriterBase):

    """Macro writer for the CMake format.

    For details on the provided methods, see MacroWriterBase, which this
    class inherits from.
    """

    def __init__(self, output):
        """Initialize a CMake macro writer.

        Arguments:
        output - File-like object (probably an io.TextIOWrapper), which
                 will be written to.
        """
        super(CMakeMacroWriter, self).__init__(output)
        # This counter is for avoiding name conflicts in temporary
        # variables used for shell commands.
        self._var_num = 0

    def environment_variable_string(self, name):
        """Return an environment variable reference.

        >>> import io
        >>> s = io.StringIO()
        >>> CMakeMacroWriter(s).environment_variable_string("foo")
        '$ENV{foo}'
        """
        return "$ENV{" + name + "}"

    def shell_command_strings(self, command):
        # pylint: disable=line-too-long
        """Return strings used to get the output of a shell command.

        >>> import io
        >>> s = io.StringIO()
        >>> set_up, inline, tear_down = CMakeMacroWriter(s).shell_command_strings("echo bar")
        >>> set_up
        'execute_process(COMMAND echo bar OUTPUT_VARIABLE CIME_TEMP_SHELL0 OUTPUT_STRIP_TRAILING_WHITESPACE)'
        >>> inline
        '${CIME_TEMP_SHELL0}'
        >>> tear_down
        'unset(CIME_TEMP_SHELL0)'
        """
        # pylint: enable=line-too-long
        # Create a unique variable name, then increment variable number
        # counter so that we get a different value next time.
        var_name = "CIME_TEMP_SHELL" + str(self._var_num)
        self._var_num += 1
        set_up = "execute_process(COMMAND " + command + \
                 " OUTPUT_VARIABLE " + var_name + \
                 " OUTPUT_STRIP_TRAILING_WHITESPACE)"
        tear_down = "unset(" + var_name + ")"
        return (set_up, "${" + var_name + "}", tear_down)

    def variable_string(self, name):
        """Return a string to refer to a variable with the given name.

        >>> import io
        >>> s = io.StringIO()
        >>> CMakeMacroWriter(s).variable_string("foo")
        '${CIME_foo}'
        """
        return "${CIME_" + name + "}"

    def set_variable(self, name, value):
        """Write out a statement setting a variable to some value.

        >>> import io
        >>> s = io.StringIO()
        >>> CMakeMacroWriter(s).set_variable("foo", "bar")
        >>> s.getvalue()
        u'set(CIME_foo "bar")\\n'
        """
        self.write_line("set(CIME_" + name + ' "' + value + '")')

    def start_ifeq(self, left, right):
        """Write out a statement to start a conditional block.

        >>> import io
        >>> s = io.StringIO()
        >>> CMakeMacroWriter(s).start_ifeq("foo", "bar")
        >>> s.getvalue()
        u'if("foo" STREQUAL "bar")\\n'
        """
        self.write_line('if("' + left + '" STREQUAL "' + right + '")')
        self.indent_right()

    def end_ifeq(self):
        """Write out a statement to end a block started with start_ifeq.

        >>> import io
        >>> s = io.StringIO()
        >>> writer = CMakeMacroWriter(s)
        >>> writer.start_ifeq("foo", "bar")
        >>> writer.set_variable("foo2", "bar2")
        >>> writer.end_ifeq()
        >>> s.getvalue()
        u'if("foo" STREQUAL "bar")\\n  set(CIME_foo2 "bar2")\\nendif()\\n'
        """
        self.indent_left()
        self.write_line("endif()")


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
            expect(not any(self.settings[i].is_ambiguous_with(other)
                           for other in self.settings[i+1:]),
                   "Variable "+self.name+" is set ambiguously in "
                   "config_build.xml. Check this file for multiple "
                   "settings that could apply to this machine.")

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
        for setting in settings:
            if setting.conditions:
                condition = setting.conditions.keys()[0]
        if condition is None:
            # If there are no conditions, we have reached a leaf.
            # We combine whatever settings are left; there should be at
            # most one non-appending setting, or an arbitrary number of
            # appending settings.
            self._is_leaf = True
            self._name = name
            self._values = []
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
                self._values.append(setting.value)
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

    def write_out(self, writer):
        """Write tree to file.

        The writer argument is an object inheriting from MacroWriterBase.
        This function first writes out all the initial settings with
        appropriate conditionals, then the appending settings.
        """
        if self._is_leaf:
            for line in self._set_up:
                writer.write_line(line)
            for value in self._values:
                if self._do_append:
                    writer.append_variable(self._name, value)
                else:
                    writer.set_variable(self._name, value)
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


class CompilerBlock(object):

    """Data used to translate a single <compiler> element.

    This is used during write_macros to traverse the XML and create a list
    of settings specified in the element.

    Public methods:
    add_settings_to_lists
    matches_machine
    """

    def __init__(self, writer, compiler_elem):
        """Construct a CompilerBlock.

        Arguments:
        writer - The Makefile/CMake writer object.
        compiler_elem - An xml.ElementTree.Element corresponding to this
                        <compiler> element.
        """
        self._writer = writer
        self._compiler_elem = compiler_elem
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
        # Create the setting object and add it to one of our
        # PossibleValues in value_lists.
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

    def matches_machine(self, os_, machine):
        """Check whether this block matches a machine/os.

        This also sets the specificity of the block, so this must be called
        before add_settings_to_lists if machine-specific output is needed.

        Arguments:
        os_ - Operating system to match.
        machine -
        """
        self._specificity = 0
        if "MACH" in self._compiler_elem.keys():
            if machine == self._compiler_elem.get("MACH"):
                self._specificity += 2
            else:
                return False
        if "OS" in self._compiler_elem.keys():
            if os_ == self._compiler_elem.get("OS"):
                self._specificity += 1
            else:
                return False
        # We allow the compiler to be changed after Macros generation,
        # so it doesn't figure into the decision of which items get
        # written to Macros.
        return True


class MacroMaker(object):

    """Class to convert config_build.xml input into a macros file.

    Public attributes:
    os - Operating system used in config_build lookup and ranking.
    machine - Machine used in config_build lookup.
    flag_vars - A set of all variables in config_build that contain "flag-
                like" data (i.e. a space-separated list of arguments).

    Public methods:
    write_macros
    """

    def __init__(self, os_, machine, schema_path=None):
        """Construct a MacroMaker given machine-specific information.

        In the process some information about possible variables is read in
        from the schema file.

        Arguments:
        os_ - Name of a machine's operating system.
        machine - Name of a machine.
        schema_path (optional) - Path to config_build.xsd within CIME.

        >>> "CFLAGS" in MacroMaker('FakeOS', 'MyMach').flag_vars
        True
        >>> "MPICC" in MacroMaker('FakeOS', 'MyMach').flag_vars
        False
        """
        self.os = os_
        self.machine = machine

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
        # pylint: disable=redefined-variable-type
        if build_system == "Makefile":
            writer = MakeMacroWriter(output)
        elif build_system == "CMake":
            writer = CMakeMacroWriter(output)
        else:
            expect(False,
                   "Unrecognized build system provided to write_macros: " +
                   build_system)
        # pylint: enable=redefined-variable-type

        # Start processing the file.
        tree = ET.parse(xml_file)
        value_lists = dict()
        for compiler_elem in tree.findall("compiler"):
            block = CompilerBlock(writer, compiler_elem)
            # If this block matches machine settings, use it.
            if block.matches_machine(self.os, self.machine):
                block.add_settings_to_lists(self.flag_vars, value_lists)

        # Now that we've scanned through the input, output the variable
        # settings.
        vars_written = set()
        while value_lists:
            made_progress = False
            # Note: We can't just iterate over value_lists because we
            # are removing entries from the dictionary as we go. Hence the
            # use of the keys method.
            for var_name in value_lists.keys():
                # If depends is a subset of vars_written, then all
                # dependencies have been handled.
                if value_lists[var_name].depends <= vars_written:
                    # Note that we're writing this variable.
                    vars_written.add(var_name)
                    made_progress = True
                    # Make the conditional trees and write them out.
                    normal_tree, append_tree = \
                                    value_lists[var_name].to_cond_trees()
                    if normal_tree is not None:
                        normal_tree.write_out(writer)
                    if append_tree is not None:
                        append_tree.write_out(writer)
                    # Remove this variable from the list and continue.
                    del value_lists[var_name]
            # If we loop through all the variables and don't write any
            # of them, then there is probably a problem with the
            # depends graph.
            expect(made_progress,
                   "The config_build XML has bad <var> references. "
                   "Check for circular references or variables that "
                   "are in a <var> tag but not actually defined.")
