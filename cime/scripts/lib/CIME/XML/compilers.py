"""
Interface to the config_compilers.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.XML.compilerblock import CompilerBlock
from CIME.BuildTools.makemacroswriter import MakeMacroWriter
from CIME.BuildTools.cmakemacroswriter import CMakeMacroWriter
from CIME.BuildTools.macroconditiontree import merge_optional_trees
import six

logger = logging.getLogger(__name__)

class Compilers(GenericXML):

    def __init__(self, machobj, infile=None, compiler=None, mpilib=None, files=None, version=None):
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
        if version is not None:
            # this is used in scripts_regression_tests to force version 2, it should not be used otherwise
            self._version = version
        else:
            self._version = self.get_version()

        self.machine  = machobj.get_machine_name()
        self.os = machobj.get_value("OS")
        if compiler is None:
            compiler = machobj.get_default_compiler()
        self.compiler       = compiler

        if mpilib is None:
            if compiler is None:
                mpilib = machobj.get_default_MPIlib()
            else:
                mpilib = machobj.get_default_MPIlib(attributes={'compiler':compiler})
        self.mpilib = mpilib

        self.compiler_nodes = None # Listed from last to first
        #Append the contents of $HOME/.cime/config_compilers.xml if it exists
        #This could cause problems if node matchs are repeated when only one is expected
        infile = os.path.join(os.environ.get("HOME"),".cime","config_compilers.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)

        if self.compiler is not None:
            self.set_compiler(compiler)

        if self._version > 1.0:
            schema_db = GenericXML(infile=schema)
            compiler_vars = schema_db.get_child("{http://www.w3.org/2001/XMLSchema}group", attributes={"name":"compilerVars"})
            choice  = schema_db.get_child(name="{http://www.w3.org/2001/XMLSchema}choice", root=compiler_vars)
            self.flag_vars = set(schema_db.get(elem, "name") for elem in schema_db.get_children(root=choice, attributes={"type":"flagsVar"}))

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
            result = self.get_optional_child(name=nodename, attributes=attributes, root=compiler_node)
            if result is not None:
                return result

        return None

    def _is_compatible(self, compiler_node, compiler, machine, os_, mpilib):
        for xmlid, value in [ ("COMPILER", compiler), ("MACH", machine), ("OS", os_), ("MPILIB", mpilib) ]:
            if value is not None and self.has(compiler_node, xmlid) and value != self.get(compiler_node, xmlid):
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
            nodes = self.get_children(name="compiler")
            for node in nodes:
                if self._is_compatible(node, compiler, machine, os_, mpilib):
                    self.compiler_nodes.append(node)

            self.compiler_nodes.reverse()

            self.compiler = compiler
            self.machine  = machine
            self.os       = os_
            self.mpilib   = mpilib

    #pylint: disable=arguments-differ
    def get_value(self, name, attribute=None, resolved=True, subgroup=None):
        """
        Get Value of fields in the config_compilers.xml file
        """
        expect(self.compiler_nodes is not None, "Compiler object has no compiler defined")
        expect(subgroup is None, "This class does not support subgroups")
        value = None

        node = self.get_optional_compiler_node(name, attributes=attribute)
        if node is not None:
            value = self.text(node)

        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

        return value

    def write_macros_file(self, macros_file="Macros.make", output_format="make", xml=None):
        if self._version <= 1.0:
            expect(False, "config_compilers.xml version '{}' is no longer supported".format(self._version))
        else:
            if output_format == "make":
                format_ = "Makefile"
            elif output_format == "cmake":
                format_ = "CMake"
            else:
                format_ = output_format

            if isinstance(macros_file, six.string_types):
                with open(macros_file, "w") as macros:
                    self._write_macros_file(format_, macros)
            else:
                self._write_macros_file(format_, macros_file, xml)

    def _write_macros_file(self, build_system, output, xml=None):
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
        node_list = []
        if xml is None:
            node_list = self.get_children(name="compiler")
        else:
            gen_xml = GenericXML()
            gen_xml.read_fd(xml)
            node_list = gen_xml.get_children(name="compiler")

        for compiler_elem in node_list:
            block = CompilerBlock(writer, compiler_elem, self._machobj, self)
            # If this block matches machine settings, use it.
            if block.matches_machine():
                block.add_settings_to_lists(self.flag_vars, value_lists)
        # Now that we've scanned through the input, output the variable
        # settings.
        vars_written = set()
        while value_lists:
            # Variables that are ready to be written.
            ready_variables = [
                var_name for var_name in value_lists
                if value_lists[var_name].dependencies() <= vars_written
            ]
            expect(len(ready_variables) > 0,
                   "The file {} has bad $VAR references. "
                   "Check for circular references or variables that "
                   "are used in a $VAR but not actually defined.".format(self.filename))
            big_normal_trees = {}
            big_append_tree = None
            for var_name in ready_variables:
                # Note that we're writing this variable.
                vars_written.add(var_name)
                # Make the conditional trees and write them out.
                normal_trees, append_tree = \
                    value_lists[var_name].to_cond_trees()
                for spec in normal_trees:
                    if spec in big_normal_trees:
                        big_normal_trees[spec] = merge_optional_trees(normal_trees[spec],
                                                                      big_normal_trees[spec])
                    else:
                        big_normal_trees[spec] = normal_trees[spec]
                big_append_tree = merge_optional_trees(append_tree,
                                                        big_append_tree)
                # Remove this variable from the list of variables to handle
                # next iteration.
                del value_lists[var_name]
            specificities = sorted(list(big_normal_trees.keys()))
            for spec in specificities:
                big_normal_trees[spec].write_out(writer)
            if big_append_tree is not None:
                big_append_tree.write_out(writer)
