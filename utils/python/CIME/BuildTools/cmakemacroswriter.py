"""Classes used to write build system files.

The classes here are used to write out settings for use by Makefile and CMake
build systems. The two relevant classes are CMakeMacroWriter and
MakeMacroWriter, which encapsulate the information necessary to write CMake and
Makefile formatted text, respectively. See the docstrings for those classes for
more.
"""

# This is not the most useful check.
# pylint: disable=invalid-name

from CIME.BuildTools.macrowriterbase import MacroWriterBase
from CIME.XML.standard_module_setup import *
logger = logging.getLogger(__name__)


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
        '${foo}'
        """
        return "${" + name + "}"

    def set_variable(self, name, value):
        """Write out a statement setting a variable to some value.

        >>> import io
        >>> s = io.StringIO()
        >>> CMakeMacroWriter(s).set_variable("foo", "bar")
        >>> s.getvalue()
        u'set(foo "bar")\\n'
        """
        print "name %s value %s"%(name,value)
        self.write_line("set(" + name + ' "' + value + '")')

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
        u'if("foo" STREQUAL "bar")\\n  set(foo2 "bar2")\\nendif()\\n'
        """
        self.indent_left()
        self.write_line("endif()")
