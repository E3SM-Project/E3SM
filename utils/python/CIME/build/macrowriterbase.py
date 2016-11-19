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

# None class based method for version 1.0

def write_macros_file_v1(macros, macros_file="Macros", output_format="make"):
    """
    Parse the config_compiler.xml file into a Macros file for the
    given machine and compiler.
    """
    # A few things can be used from environ if not in XML
    for item in ["MPI_PATH", "NETCDF_PATH"]:
	if not item in macros and item in os.environ:
	    logger.warn("Setting %s from Environment" % item)
	    macros[item] = os.environ[item]

    with open(macros_file, "w") as fd:
	fd.write(
"""#
# COMPILER=%s
# OS=%s
# MACH=%s
""" % (self.compiler, self.os, self.machine)
)
	if output_format == "make":
	    fd.write("#\n# Makefile Macros generated from %s \n#\n" % self.filename)

	    # print the settings out to the Macros file
	    for key, value in sorted(macros.iteritems()):
		if key == "_COND_":
		    pass
		elif key.startswith("ADD_"):
		    fd.write("%s+=%s\n\n" % (key[4:], value))
		else:
		    fd.write("%s:=%s\n\n" % (key, value))

	elif output_format == "cmake":
	    fd.write(
'''#
# cmake Macros generated from $compiler_file
#
include(Compilers)
set(CMAKE_C_FLAGS_RELEASE "" CACHE STRING "Flags used by c compiler." FORCE)
set(CMAKE_C_FLAGS_DEBUG "" CACHE STRING "Flags used by c compiler." FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "" CACHE STRING "Flags used by Fortran compiler." FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "" CACHE STRING "Flags used by Fortran compiler." FORCE)
set(all_build_types "None Debug Release RelWithDebInfo MinSizeRel")
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build, options are: ${all_build_types}." FORCE)
''')

	    # print the settings out to the Macros file, do it in
	    # two passes so that path values appear first in the
	    # file.
	    for key, value in sorted(macros.iteritems()):
		if key == "_COND_":
		    pass
		else:
		    value = value.replace("(", "{").replace(")", "}")
		    if key.endswith("_PATH"):
			fd.write("set(%s %s)\n" % (key, value))
			fd.write("list(APPEND CMAKE_PREFIX_PATH %s)\n\n" % value)

	    for key, value in sorted(macros.iteritems()):
		if key == "_COND_":
		    pass
		else:
		    value = value.replace("(", "{").replace(")", "}")
		    if "CFLAGS" in key:
			fd.write("add_flags(CMAKE_C_FLAGS %s)\n\n" % value)
		    elif "FFLAGS" in key:
			fd.write("add_flags(CMAKE_Fortran_FLAGS %s)\n\n" % value)
		    elif "CPPDEFS" in key:
			fd.write("list(APPEND COMPILE_DEFINITIONS %s)\n\n" % value)
		    elif "SLIBS" in key or "LDFLAGS" in key:
			fd.write("add_flags(CMAKE_EXE_LINKER_FLAGS %s)\n\n" % value)

	# Recursively print the conditionals, combining tests to avoid repetition
	_parse_hash(macros["_COND_"], fd, 0, output_format)


def _parse_hash(macros, fd, depth, output_format, cmakedebug=""):
    width = 2 * depth
    for key, value in macros.iteritems():
	if type(value) is dict:
	    if output_format == "make" or "DEBUG" in key:
		for key2, value2 in value.iteritems():
		    if output_format == "make":
			fd.write("%sifeq ($(%s), %s) \n" % (" " * width, key, key2))

		    _parse_hash(value2, fd, depth + 1, output_format, key2)
	else:
	    if output_format == "make":
		if key.startswith("ADD_"):
		    fd.write("%s %s += %s\n" % (" " * width, key[4:], value))
		else:
		    fd.write("%s %s += %s\n" % (" " * width, key, value))

	    else:
		value = value.replace("(", "{").replace(")", "}")
		release = "DEBUG" if "TRUE" in cmakedebug else "RELEASE"
		if "CFLAGS" in key:
		    fd.write("add_flags(CMAKE_C_FLAGS_%s %s)\n\n" % (release, value))
		elif "FFLAGS" in key:
		    fd.write("add_flags(CMAKE_Fortran_FLAGS_%s %s)\n\n" % (release, value))
		elif "CPPDEF" in key:
		    fd.write("add_config_definitions(%s %s)\n\n" % (release, value))
		elif "SLIBS" in key or "LDFLAGS" in key:
		    fd.write("add_flags(CMAKE_EXE_LINKER_FLAGS_%s %s)\n\n" % (release, value))

    width -= 2
    if output_format == "make" and depth > 0:
	fd.write("%sendif\n\n" % (" " * width))

def _add_to_macros(node, macros):
    for child in node:
	name = child.tag
	attrib = child.attrib
	value = child.text

	if not attrib:
	    if name.startswith("ADD_"):
		basename = name[4:]
		if basename in macros:
		    macros[basename] = "%s %s" % (macros[basename], value)
		elif name in macros:
		    macros[name] = "%s %s" % (macros[name], value)
		else:
		    macros[name] = value
	    else:
		macros[name] = value

	else:
	    cond_macros = macros["_COND_"]
	    for key, value2 in attrib.iteritems():
		if key not in cond_macros:
		    cond_macros[key] = {}
		if value2 not in cond_macros[key]:
		    cond_macros[key][value2] = {}
		cond_macros = cond_macros[key][value2]

	    cond_macros[name] = value
