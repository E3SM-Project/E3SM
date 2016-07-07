"""Classes used to write build system files.

The classes here are used to write out settings for use by Makefile and CMake
build systems. The two relevant classes are CMakeMacroWriter and
MakeMacroWriter, which encapsulate the information necessary to write CMake and
Makefile formatted text, respectively. See the docstrings for those classes for
more.
"""

# This is not the most useful check.
# pylint: disable=invalid-name

from abc import ABCMeta, abstractmethod

__all__ = ["CMakeMacroWriter", "MakeMacroWriter"]

class MacroWriterBase(object):

    """Abstract base class for macro file writers.

    The methods here come in three flavors:
    1. indent_left/indent_right change the level of indent used internally by
       the class.
    2. The various methods ending in "_string" return strings relevant to the
       build system.
    3. The other methods write information to the file handle associated with
       an individual writer instance.

    Public attributes:
    indent_increment - Number of spaces to indent if blocks (does not apply
                       to format-specific indentation, e.g. cases where
                       Makefiles must use tabs).
    output - File-like object that output is written to.

    Public methods:
    indent_string
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

    def indent_string(self):
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
        self.output.write(unicode(self.indent_string() + line + "\n"))

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
