"""Functions to set up machine-specific settings.

Currently, this module is not under unit test, because it interacts too
directly with the OS of the machine it is on. Any changes to address this
would be welcome (this would likely require the use of mock modules,
and/or a specific set of machines to test on).

Public classes:
MachineCompilerSettings - Set up machine/compiler specific environment.
"""

# Python 3 compatible printing in Python 2.
from __future__ import print_function

from itertools import chain
from os import environ
import platform
import re
import subprocess
from xml.etree.ElementTree import ElementTree

import environment as ENV
from printer import ScriptPrinter
from xml_utils import best_match, all_matches

__all__ = ("MachineCompilerSettings")


def get_machine_name():
    """Returns a canonical version of the machine name."""
    name = platform.node()
    # Strip everything after the first ".", and whitespace.
    name = name.split(".")[0].strip()
    if re.match("^yslogin[0-9]+", name):
        name = "yellowstone"
    elif re.match("^miralac[0-9]+", name):
        name = "mira"
    return name

def load_machine_env(compiler):
    """Add machine environment variables not in config_compilers.xml

    Besides simply setting variables, this may also load some modules.
    """

    mach = get_machine_name()

    if mach == "yellowstone":
        mod = ENV.ModuleInterface()
        mod.python_init("/glade/apps/opt/lmod/lmod/init/env_modules_python.py")
        mod.purge()
        mod.load("ncarenv/1.0")
        mod.load("ncarbinlibs/1.0")
        if compiler == "intel":
            mod.load("intel/13.1.2")
        elif compiler == "pgi":
            mod.load("pgi/13.9")
        mod.load("ncarcompilers/1.0")
        mod.load("cmake/2.8.10.2")


class MachineCompilerSettings(object):

    """Encapsulate machine settings and set environment from them.

    Public methods:
    __init__ - Discover information about local machine.
    compiler_xml_to_env - Set a specific variable using config_compilers.
    set_compiler_env - Set up machine-specific environment.
    write_cmake_macros - Create CMake "Macros" file.
    """

    def __init__(self, compiler, compiler_xml_path, use_mpi=False):
        """Discover information about the machine and compiler.

        Arguments:
        compiler - String naming the compiler vendor.
        compiler_xml_path - Path to config_compilers.xml file.
        use_mpi - Boolean option to turn on use of MPI compiler wrappers.
        """
        self.machine_dict = {
            "COMPILER": compiler,
            "MACH": get_machine_name(),
            "OS": platform.system(),
            }
        self.compiler_xml_tree = ElementTree()
        self.compiler_xml_tree.parse(compiler_xml_path)
        self.use_mpi = use_mpi

    def compiler_xml_to_env(self, xml_path, var_name):
        """Look up a config_compilers entry and set a variable from it.

        Arguments:
        xml_path - Path within the xml file (e.g. "compiler/SFC").
        var_name - Name of environment variable (e.g. "FC").
        """
        match = best_match(self.compiler_xml_tree, xml_path,
                           self.machine_dict)
        assert match is not None, "Could not determine "+var_name+ \
            " from compiler/machines xml file."
        environ[var_name] = match.text.strip()

    def set_compiler_env(self):
        """Set up the environment for this machine and compiler."""
        if self.use_mpi:
            self.compiler_xml_to_env("compiler/MPIFC", "FC")
            self.compiler_xml_to_env("compiler/MPICC", "CC")
        else:
            self.compiler_xml_to_env("compiler/SFC", "FC")
            self.compiler_xml_to_env("compiler/SCC", "CC")
        load_machine_env(self.machine_dict["COMPILER"])

    def write_cmake_macros(self, macros_file):
        """Write CMake macros file using config_compilers.xml

        Arguments:
        macros_file - File object to write to.
        """

        # Print header to file.
        macros_printer = ScriptPrinter(macros_file)
        header_lines = [
            "CESM build flags for:",
            "  Compiler = "+self.machine_dict["COMPILER"],
            "  Machine = "+self.machine_dict["MACH"],
            "  OS = "+self.machine_dict["OS"],
            ]
        for line in header_lines:
            macros_printer.comment(line)

        # pFUnit location if it exists.
        match = best_match(self.compiler_xml_tree, "compiler/PFUNIT_PATH",
                           self.machine_dict)
        if match is not None:
            macros_printer.print_header("pFUnit location.")
            macros_printer.print(
                "list(APPEND CMAKE_PREFIX_PATH "+match.text+")"
                )

        # Normal and debug dictionaries for looking things up in
        # config_compilers.
        normal_dict = self.machine_dict.copy()
        normal_dict["DEBUG"] = "FALSE"

        debug_dict = self.machine_dict.copy()
        debug_dict["DEBUG"] = "TRUE"

        def add_formatted_flags(flags_name, format):
            """Print CMake flags using macros_printer.

            Arguments:
            flags_name - Name to search for in config_compilers.
            format - Function that takes a build type and flag match, and
            returns the string to print out.
            """

            paths = ["compiler/"+flags_name, "compiler/ADD_"+flags_name]

            # This creates an iterable over elements in config_compilers
            # that match in non-debug mode.
            normal_matches = chain.from_iterable(
                all_matches(self.compiler_xml_tree, path, normal_dict)
                for path in paths
                )
            for match in normal_matches:
                macros_printer.print(format("CESM", match.text))

            # Now the same for debug mode.
            debug_matches = chain.from_iterable(
                all_matches(self.compiler_xml_tree, path, debug_dict)
                for path in paths
                )
            for match in debug_matches:
                macros_printer.print(format("CESM_DEBUG", match.text))

        # Below, we just use a bunch of lambda functions to describe how
        # the build type and a matching element (e.g. an FFLAGS entry) are
        # turned into a CMake function call.

        macros_printer.print_header("CPP definitions.")
        add_formatted_flags(
            "CPPDEFS",
            lambda b, m: "add_config_definitions("+b+" "+m+")"
            )
        def format_contiguous(build_type, match):
            comma = "," if self.machine_dict["COMPILER"] != "ibm" else "\\\,"
            contig_def = "contiguous"+comma if match == "TRUE" else ""
            return "add_config_definitions("+build_type+\
                " -DUSE_CONTIGUOUS="+contig_def+")"
        add_formatted_flags(
            "HAS_F2008_CONTIGUOUS",
            format_contiguous
            )

        macros_printer.print_header("Fortran flags.")
        add_formatted_flags(
            "FFLAGS",
            lambda b, m: "add_flags(CMAKE_Fortran_FLAGS_"+b+" "+m+")"
            )

        macros_printer.print_header("C flags.")
        add_formatted_flags(
            "CFLAGS",
            lambda b, m: "add_flags(CMAKE_C_FLAGS_"+b+" "+m+")"
            )

        macros_printer.print_header("Linker flags.")
        add_formatted_flags(
            "LDFLAGS",
            lambda b, m: "add_flags(CMAKE_EXE_LINKER_FLAGS_"+b+" "+m+")"
            )
