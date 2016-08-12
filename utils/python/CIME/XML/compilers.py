"""
Interface to the config_compilers.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files

logger = logging.getLogger(__name__)

class Compilers(GenericXML):

    def __init__(self, compiler=None, machine=None, os_= None, mpilib=None, infile=None, files=None):
        """
        initialize an object
        """
        if infile is None:
            if files is None:
                files = Files()
            infile = files.get_value("COMPILERS_SPEC_FILE")

        GenericXML.__init__(self, infile)

        self.machine        = machine
        self.os             = os_
        self.mpilib         = mpilib
        self.compiler_nodes = None # Listed from last to first
        self.compiler       = compiler
        #Append the contents of $HOME/.cime/config_compilers.xml if it exists
        #This could cause problems if node matchs are repeated when only one is expected
        infile = os.path.join(os.environ.get("HOME"),".cime","config_compilers.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)

        if self.compiler is not None:
            self.set_compiler(compiler)


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

        >>> machobj = Compilers(machine="melvin")
        >>> machobj.set_compiler("gnu")
        >>> machobj.get_compiler()
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
        """
        Parse the config_compiler.xml file into a Macros file for the
        given machine and compiler.
        """

        # Parse the xml settings into the $macros hash structure
        # put conditional settings in the _COND_ portion of the hash
        # and handle them seperately
        macros = {"_COND_" : {}}

        # Do worst matches first
        for compiler_node in reversed(self.compiler_nodes):
            _add_to_macros(compiler_node, macros)

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
