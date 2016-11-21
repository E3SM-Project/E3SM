"""
Interface to the config_compilers.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.XML.compilerblock import CompilerBlock
from CIME.BuildTools.macrowriterbase import write_macros_file_v1
from CIME.BuildTools.makemacroswriter import MakeMacroWriter
from CIME.BuildTools.cmakemacroswriter import CMakeMacroWriter
from CIME.BuildTools.macroconditiontree import merge_optional_trees

logger = logging.getLogger(__name__)

class Compilers(GenericXML):

    def __init__(self, machobj, infile=None, compiler=None, mpilib=None, files=None):
        """
        initialize an object
        """
        if infile is None:
            if files is None:
                files = Files()
            infile = files.get_value("COMPILERS_SPEC_FILE")
            schema = files.get_schema("COMPILERS_SPEC_FILE")

        GenericXML.__init__(self, infile, schema)
        self._machobj = machobj
        self._version = self.get_version()

        self.machine  = machobj.get_machine_name()
        self.os = machobj.get_value("OS")
        if mpilib is None:
            mpilib = machobj.get_default_MPIlib()
        self.mpilib = mpilib
        if compiler is None:
            compiler = machobj.get_default_compiler()
        self.compiler       = compiler

        self.compiler_nodes = None # Listed from last to first
        #Append the contents of $HOME/.cime/config_compilers.xml if it exists
        #This could cause problems if node matchs are repeated when only one is expected
        infile = os.path.join(os.environ.get("HOME"),".cime","config_compilers.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)

        if self.compiler is not None:
            self.set_compiler(compiler)

        if self._version != "1.0":
            # Run an XPath query to extract the list of flag variable names.
            ns = {"xs": "http://www.w3.org/2001/XMLSchema"}
            flag_xpath = ".//xs:group[@name='compilerVars']/xs:choice/xs:element[@type='flagsVar']"
            flag_elems = ET.parse(schema).getroot().findall(flag_xpath, ns)
            self.flag_vars = set(elem.get('name') for elem in flag_elems)



    def get_compiler(self):
        """
        Return the name of the compiler
        """
        return self.compiler

    def get_optional_compiler_node(self, nodename, attributes=None):
        """
        Return data on a node for a compiler
        """
        expect(self.compiler_nodes is not None, "Compiler not set, use parent get_node?")
        for compiler_node in self.compiler_nodes:
            result = self.get_optional_node(nodename, attributes, root=compiler_node)
            if result is not None:
                return result

        return None

    def _is_compatible(self, compiler_node, compiler, machine, os_, mpilib):
        for xmlid, value in [ ("COMPILER", compiler), ("MACH", machine), ("OS", os_), ("MPILIB", mpilib) ]:
            if value is not None and xmlid in compiler_node.attrib and value != compiler_node.get(xmlid):
                return False

        return True

    def set_compiler(self, compiler, machine=None, os_=None, mpilib=None):
        """
        Sets the compiler block in the Compilers object
        >>> from CIME.XML.machines import Machines
        >>> compobj = Compilers(Machines(machine="melvin"))
        >>> compobj.set_compiler("gnu")
        >>> compobj.get_compiler()
        'gnu'
        """
        machine = machine if machine else self.machine
        os_     = os_ if os_ else self.os
        mpilib  = mpilib if mpilib else self.mpilib

        if self.compiler != compiler or self.machine != machine or self.os != os_ or self.mpilib != mpilib or self.compiler_nodes is None:
            self.compiler_nodes = []
            nodes = self.get_nodes("compiler")
            for node in nodes:
                if self._is_compatible(node, compiler, machine, os_, mpilib):
                    self.compiler_nodes.append(node)

            self.compiler_nodes.reverse()

            self.compiler = compiler
            self.machine  = machine
            self.os       = os_
            self.mpilib   = mpilib

    def get_value(self, name, attribute=None, resolved=True, subgroup=None):
        """
        Get Value of fields in the config_compilers.xml file
        """
        expect(self.compiler_nodes is not None, "Compiler object has no compiler defined")
        expect(subgroup is None, "This class does not support subgroups")
        value = None

        node = self.get_optional_compiler_node(name)
        if node is not None:
            value = node.text

        if value is None:
            # if all else fails
            value = GenericXML.get_value(self, name)

        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

        return value

    def write_macros_file(self, macros_file="Macros", output_format="make"):
        if self._version == "1.0":
            # Parse the xml settings into the $macros hash structure
            # put conditional settings in the _COND_ portion of the hash
            # and handle them seperately
            macros = {"_COND_" : {}}

            # Do worst matches first
            for compiler_node in reversed(self.compiler_nodes):
                _add_to_macros(compiler_node, macros)
            write_macros_file_v1(macros, self.compiler, self.os,
                                        self.machine, macros_file,
                                        output_format)
        else:
            if output_format == "make":
                format = "Makefile"
            elif output_format == "cmake":
                format = "CMake"
            with open(macros_file, "w") as macros:
                self._write_macros_file_v2(format, macros)

    def _write_macros_file_v2(self, build_system, output):
        """Write a Macros file for this machine.

        Arguments:
        build_system - Format of the file to be written. Currently the only
                       valid values are "Makefile" and "CMake".
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
        for compiler_elem in self.get_nodes("compiler"):
            block = CompilerBlock(writer, compiler_elem, self._machobj)
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
                big_normal_tree = merge_optional_trees(normal_tree,
                                                        big_normal_tree)
                big_append_tree = merge_optional_trees(append_tree,
                                                        big_append_tree)
                # Remove this variable from the list of variables to handle
                # next iteration.
                del value_lists[var_name]
            if big_normal_tree is not None:
                big_normal_tree.write_out(writer)
            if big_append_tree is not None:
                big_append_tree.write_out(writer)



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
