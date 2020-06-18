"""Classes used to write build system files.

The classes here are used to write out settings for use by Makefile and CMake
build systems. The two relevant classes are CMakeMacroWriter and
MakeMacroWriter, which encapsulate the information necessary to write CMake and
Makefile formatted text, respectively. See the docstrings for those classes for
more.
"""

from CIME.BuildTools.macrowriterbase import MacroWriterBase
from CIME.XML.standard_module_setup import *
#logger = logging.getLogger(__name__)

# This is not the most useful check.
# pylint: disable=invalid-name

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
        >>> str(s.getvalue())
        'foo := bar\\n'
        """
        # Note that ":=" is used so that we can control the behavior for
        # both Makefile and CMake variables similarly.
        self.write_line(name + " := " + value)

    def start_ifeq(self, left, right):
        """Write out a statement to start a conditional block.

        >>> import io
        >>> s = io.StringIO()
        >>> MakeMacroWriter(s).start_ifeq("foo", "bar")
        >>> str(s.getvalue())
        'ifeq (foo,bar)\\n'
        """
        if right.startswith("!"):
            right = right.lstrip("!")
            not_str = "n"
        else:
            not_str = ""

        self.write_line("if{}eq ({},{})".format(not_str, left, right))
        self.indent_right()

    def end_ifeq(self):
        """Write out a statement to end a block started with start_ifeq.

        >>> import io
        >>> s = io.StringIO()
        >>> writer = MakeMacroWriter(s)
        >>> writer.start_ifeq("foo", "bar")
        >>> writer.set_variable("foo2", "bar2")
        >>> writer.end_ifeq()
        >>> str(s.getvalue())
        'ifeq (foo,bar)\\n  foo2 := bar2\\nendif\\n'
        """
        self.indent_left()
        self.write_line("endif")
