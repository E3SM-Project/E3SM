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
from CIME.XML.standard_module_setup import *
from CIME.utils import get_cime_root
from six import add_metaclass

logger = logging.getLogger(__name__)

def _get_components(value):
    """
    >>> value = '-something ${shell ${NETCDF_PATH}/bin/nf-config --flibs} -lblas -llapack'
    >>> _get_components(value)
    [(False, '-something'), (True, '${NETCDF_PATH}/bin/nf-config --flibs'), (False, '-lblas -llapack')]
    >>> value = '${shell ${NETCDF_PATH}/bin/nf-config --flibs} -lblas -llapack'
    >>> _get_components(value)
    [(True, '${NETCDF_PATH}/bin/nf-config --flibs'), (False, '-lblas -llapack')]
    >>> value = '${shell ${NETCDF_PATH}/bin/nf-config --flibs}'
    >>> _get_components(value)
    [(True, '${NETCDF_PATH}/bin/nf-config --flibs')]
    """
    value = value.strip()
    components = []
    curr_comp = ""
    idx = 0
    while idx < len(value):
        if value[idx:idx+8] == "${shell ":
            if curr_comp:
                components.append((False, curr_comp.strip()))
                curr_comp = ""

            idx += 8
            brace_cnt = 0
            done = False
            while not done:
                if value[idx] == "{":
                    brace_cnt += 1
                    curr_comp += value[idx]

                elif value[idx] == "}":
                    if brace_cnt == 0:
                        done = True
                    else:
                        brace_cnt -= 1
                        curr_comp += value[idx]

                else:
                    curr_comp += value[idx]

                idx += 1

            components.append((True, curr_comp.strip()))
            curr_comp = ""
        else:
            curr_comp += value[idx]
            idx += 1

    if curr_comp:
        components.append((False, curr_comp.strip()))

    return components

@add_metaclass(ABCMeta)
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
        self.output.write(u"{}{}\n".format(self.indent_string(), line))

    @abstractmethod
    def environment_variable_string(self, name):
        """Return an environment variable reference."""

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

    @abstractmethod
    def variable_string(self, name):
        """Return a string to refer to a variable with the given name."""

    @abstractmethod
    def set_variable(self, name, value):
        """Write out a statement setting a variable to some value."""

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

    @abstractmethod
    def end_ifeq(self):
        """Write out a statement to end a block started with start_ifeq."""
