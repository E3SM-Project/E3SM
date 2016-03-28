"""
API for setup?
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.compilers import Compilers
from CIME.case import Case

import glob

logger = logging.getLogger(__name__)

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
            for key, value in attrib.iteritems():
                if key not in cond_macros:
                    cond_macros[key] = {}
                if value not in cond_macros[key]:
                    cond_macros[key][value] = {}
                cond_macros = cond_macros[key][value]

            cond_macros[name] = value

def _parse_hash(macros, fd, depth, output_format, cmakedebug=""):
    width = 2 * depth
    for key, value in macros.iteritems():
        if type(value) is dict:
            if output_format == "make" or "DEBUG" in key:
                for key2, value2 in value.iteritems():
                    if output_format == "make":
                        fd.write("%sifeq ($(%s), %s) \n" % (" " * width, key, key2))

                    _parse_hash(value2, depth + 1, output_format, key2)
        else:
            if output_format == "make":
                if key.startswith("ADD_"):
                    fd.write("%s %s += %s\n", " " * width, key[4:], value)
                else:
                    fd.write("%s %s += %s\n", " " * width, key, value)

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

def set_compiler(compiler_file_arg=None, case=None, macros_file="Macros", output_format="make"):
    """
    Parse the config_compiler.xml file into a Macros file for the
    given machine and compiler. Search the user's ~/.cime directory
    first, then use the standard compiler file if it is not available.
    """
    case = Case() if case is None else case

    compiler = case.get_value("COMPILER")
    machine  = case.get_value("MACH")
    os_      = case.get_value("OS")
    mpilib   = case.get_value("MPILIB")

    compilers = None
    # JGF: Which should take precedence?
    for compiler_file in [os.path.join(os.environ["HOME"], ".cime", "config_compilers.xml"),
                          compiler_file_arg]:
        if compiler_file is None or os.path.exists(compiler_file):
            compilers = Compilers(compiler=compiler, machine=machine, os_=os_, mpilib=mpilib, infile=compiler_file)
            if compiler_file is None:
                compiler_file = compilers.filename
            break

    expect(compilers is not None, "File does not exist: %s" % compiler_file_arg)

    # Parse the xml settings into the $macros hash structure
    # put conditional settings in the _COND_ portion of the hash
    # and handle them seperately
    macros = {"_COND_" : {}}

    # Do parent first
    if (compilers.parent_node is not compilers.compiler_node):
        _add_to_macros(compilers.parent_node, macros)

    _add_to_macros(compilers.compiler_node, macros)

    compcpp = compiler.upper()
    macros["ADD_CPPDEFS"] += " -D%s -DCPR%s " % (os, compcpp)

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
""" % (compiler, os_, machine)
)
        if output_format == "make":
            fd.write("#\n# Makefile Macros generated from %s \n#\n" % compiler_file)

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
