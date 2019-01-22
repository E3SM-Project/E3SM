"""
Classes used to build the CIME Macros file.

The main "public" class here is Build. It is initialized with machine-specific
information, and its write_macros method is the driver for translating the
config_compilers.xml file into a Makefile or CMake-format Macros file.

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
  $VAR, $ENV{...} and $SHELL{...} and keeping track of dependencies induced by one
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

from CIME.XML.standard_module_setup import *
from CIME.BuildTools.valuesetting import ValueSetting
from CIME.BuildTools.possiblevalues import PossibleValues

logger = logging.getLogger(__name__)

class CompilerBlock(object):

    """Data used to translate a single <compiler> element.

    This is used during write_macros to traverse the XML and create a list
    of settings specified in the element.

    Public methods:
    add_settings_to_lists
    matches_machine
    """

    def __init__(self, writer, compiler_elem, machobj, db):
        """Construct a CompilerBlock.

        Arguments:
        writer - The Makefile/CMake writer object.
        compiler_elem - An xml.ElementTree.Element corresponding to this
                        <compiler> element.
        machobj - Machines object for this machine.
        """
        self._writer = writer
        self._compiler_elem = compiler_elem
        self._db            = db
        self._machobj = machobj
        # If there's no COMPILER attribute, self._compiler is None.
        self._compiler = db.get(compiler_elem, "COMPILER")
        self._specificity = 0

    def _handle_references(self, elem, set_up, tear_down, depends):
        """Expand markup used internally.

        This function is responsible for expanding $ENV{...}, $VAR, and
        $SHELL{...} syntax into Makefile/CMake syntax.

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
        output = self._db.text(elem)
        if output is None:
            output = ""

        logger.debug("Initial output={}".format(output))
        reference_re = re.compile(r'\${?(\w+)}?')
        env_ref_re   = re.compile(r'\$ENV\{(\w+)\}')
        shell_prefix = "$SHELL{"

        for m in reference_re.finditer(output):
            var_name = m.groups()[0]
            if var_name not in ("SHELL","ENV"):
                output = output.replace(m.group(), writer.variable_string(var_name))
                depends.add(var_name)

        logger.debug("preenv pass output={}".format(output))

        for m in env_ref_re.finditer(output):
            logger.debug("look for {} in env {}".format(output,writer.environment_variable_string(m.groups()[0])))
            output = output.replace(m.group(),
                                    writer.environment_variable_string(m.groups()[0]))
            logger.debug("and output {}".format(output))

        logger.debug("postenv pass output={}".format(output))

        while shell_prefix in output:
            sidx = output.index(shell_prefix)
            brace_count = 1
            idx = 0
            for idx in range(sidx + len(shell_prefix), len(output)):
                if output[idx] == "{":
                    brace_count += 1
                elif output[idx] == "}":
                    brace_count -= 1
                    if brace_count == 0:
                        break

            command = output[sidx + len(shell_prefix) : idx]
            logger.debug("execute {} in shell, command {}".format(output, command))
            new_set_up, inline, new_tear_down = \
                writer.shell_command_strings(command)
            output = output.replace(output[sidx:idx+1], inline, 1)
            if new_set_up is not None:
                set_up.append(new_set_up)
            if new_tear_down is not None:
                tear_down.append(new_tear_down)
            logger.debug("set_up {} inline {} tear_down {}".format(new_set_up,inline,new_tear_down))

        logger.debug("First pass output={}".format(output))

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
        conditions = self._db.attrib(elem)
        if self._compiler is not None:
            conditions["COMPILER"] = self._compiler
        # Deal with internal markup.
        set_up = []
        tear_down = []
        depends = set()
        value_text = self._handle_references(elem, set_up,
                                             tear_down, depends)
        # Create the setting object.
        append = self._db.name(elem) == "append"
        setting = ValueSetting(value_text, append,
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
            value_lists[name].add_setting(setting, self._specificity,depends)

    def add_settings_to_lists(self, flag_vars, value_lists):
        """Add all data in the <compiler> element to lists of settings.

        Arguments:
        flag_vars - A set of variables containing "flag-like" data.
        value_lists - A dictionary of PossibleValues, containing the lists
                      of all settings for each variable.
        """
        for elem in self._db.get_children(root=self._compiler_elem):
            # Deal with "flag"-type variables.
            if self._db.name(elem) in flag_vars:
                for child in self._db.get_children(root=elem):
                    self._add_elem_to_lists(self._db.name(elem), child, value_lists)
            else:
                self._add_elem_to_lists(self._db.name(elem), elem, value_lists)

    def matches_machine(self):
        """Check whether this block matches a machine/os.
        This also sets the specificity of the block, so this must be called
        before add_settings_to_lists if machine-specific output is needed.
        """
        self._specificity = 0
        if self._db.has(self._compiler_elem, "MACH"):
            if self._machobj.get_machine_name() == \
               self._db.get(self._compiler_elem, "MACH"):
                self._specificity += 2
            else:
                return False
        if self._db.has(self._compiler_elem, "OS"):
            if self._machobj.get_value("OS") == self._db.get(self._compiler_elem, "OS"):
                self._specificity += 1
            else:
                return False
        # Check if the compiler is valid on this machine.
        if self._compiler is not None:
            return self._machobj.is_valid_compiler(self._compiler)
        else:
            return True
