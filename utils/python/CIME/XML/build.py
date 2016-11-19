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
from CIME.XML.standard_module_setup import *

__all__ = ["Build"]

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
