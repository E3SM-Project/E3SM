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
import os.path
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
    elif re.match("^cetuslac[0-9]+", name):
        name = "mira"
    elif re.match("^caldera.*", name) or re.match("^pronghorn.*", name):
        # Use yellowstone settings for caldera/pronghorn, since the mahcines
        # files aren't set up explicitly for those machines, and they have the
        # same configuration as yellowstone.
        name = "yellowstone"
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
            mod.load("intel/15.0.1")
            mod.load("mkl/11.1.2")
        elif compiler == "pgi":
            mod.load("pgi/13.9")
        mod.load("ncarcompilers/1.0")
        mod.load("cmake/2.8.10.2")
        mod.load("netcdf/4.3.0")


class MachineCompilerSettings(object):

    """Encapsulate machine settings and set environment from them.

    Public methods:
    __init__ - Discover information about local machine.
    compiler_xml_to_env - Set a specific variable using config_compilers.
    set_compiler_env - Set up machine-specific environment.
    write_cmake_macros - Create CMake "Macros" file.
    """

    def __init__(self, compiler, compiler_xml_path, 
                 machine=None,
                 use_env_compiler=False,
                 mpilib=None, use_openmp=False):
        """Discover information about the machine and compiler.

        Arguments:
        compiler - String naming the compiler vendor.
        compiler_xml_path - Path to config_compilers.xml file.
        machine - String naming the machine (guessed if not provided)
        use_env_compiler - Force use of environment variable FC instead of
                           using the compiler in config_compilers.xml.
        mpilib - String naming the mpi library if any.
        use_openmp - Boolean option to use OpenMP compiler flags.
                     Currently only works for CESM/CESM_DEBUG builds.
        """
        if machine is None:
            machine = get_machine_name()
        self.machine_dict = {
            "COMPILER": compiler,
            "MACH": machine,
            "OS": platform.system(),
            "compile_threaded": "true" if use_openmp else "false",
            }
        if mpilib:
            self.machine_dict["MPILIB"] = mpilib
        else:
            self.machine_dict["MPILIB"] = 'mpi-serial'

        self.compiler_xml_tree = ElementTree()
        self.compiler_xml_tree.parse(compiler_xml_path)
        self.use_env_compiler = use_env_compiler
        self.mpilib = mpilib

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

    def add_path(self, name, macros_printer):
        match = best_match(self.compiler_xml_tree, "compiler/"+name+"_PATH",
                           self.machine_dict)
        if match is not None:
            macros_printer.print_header(name + " location.")
            libpath = match.text
            _make_env_re = re.compile(
                """\$\(                    # Initial "$" and brace
                (?P<name>[A-Za-z0-9_]+) # Variable name
                \)                      # Close brace""", re.X)
            libpath = _make_env_re.sub( "$ENV{\g<name>}",libpath)
            print("libpath = "+libpath)
            macros_printer.print(
                "list(APPEND CMAKE_PREFIX_PATH "+libpath+")"
                )

    def set_compiler_env(self):
        """Set up the environment for this machine and compiler."""
        if not self.use_env_compiler:
            if (self.mpilib == "mpi-serial" or self.mpilib is None):
                self.compiler_xml_to_env("compiler/SFC", "FC")
                self.compiler_xml_to_env("compiler/SCC", "CC")
            else:
                self.compiler_xml_to_env("compiler/MPIFC", "FC")
                self.compiler_xml_to_env("compiler/MPICC", "CC")
        load_machine_env(self.machine_dict["COMPILER"])

    def write_cmake_macros(self, macros_file, model="CESM"):
        """Write CMake macros file using config_compilers.xml

        Arguments:
        macros_file - File object to write to.
        """

        # Print header to file.
        macros_printer = ScriptPrinter(macros_file)
        header_lines = [model+
            " build flags for:",
            "  Compiler = "+self.machine_dict["COMPILER"],
            "  Machine = "+self.machine_dict["MACH"],
            "  OS = "+self.machine_dict["OS"],
            " MPILIB = "+self.machine_dict["MPILIB"],
            ]

        for line in header_lines:
            macros_printer.comment(line)

        macros_printer.print("include(Compilers)")
        macros_printer.print(
            "set(CMAKE_Fortran_FLAGS_"+model+" \"\" CACHE STRING \"Flags used by the Fortran compiler during builds.\" FORCE)")
        macros_printer.print(
            "set(CMAKE_Fortran_FLAGS_"+model+"_DEBUG \"\" CACHE STRING \"Flags used by the Fortran compiler during DEBUG builds.\" FORCE)")
            

        macros_printer.print(
            "set(CMAKE_C_FLAGS_"+model+" \"\" CACHE STRING \"Flags used by the C compiler during builds.\" FORCE)")
        macros_printer.print(
            "set(CMAKE_C_FLAGS_"+model+"_DEBUG \"\" CACHE STRING \"Flags used by the C compiler during DEBUG builds.\" FORCE)")
            
        macros_printer.print(
            "mark_as_advanced(CMAKE_Fortran_FLAGS_"+model+" CMAKE_Fortran_FLAGS_"
            +model+"_DEBUG)")
        macros_printer.print(
            "mark_as_advanced(CMAKE_C_FLAGS_"+model+" CMAKE_C_FLAGS_"+model+"_DEBUG)")
        macros_printer.print(
            "set(all_build_types \"None Debug Release RelWithDebInfo MinSizeRel "+model+
            " "+model+"_DEBUG\")")
        macros_printer.print(
            "set(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\" CACHE STRING \"Choose the type of build, options are: ${all_build_types}.\" FORCE)")
            



        # pFUnit location if it exists.
        self.add_path("PFUNIT",macros_printer)
        # NETCDF location if it exists.
        self.add_path("NETCDF",macros_printer)
        # PNETCDF location if it exists.
        self.add_path("PNETCDF", macros_printer)
        # HDF5 location if it exists.
        self.add_path("HDF5", macros_printer)
        # MPI location if it exists.
        self.add_path("MPI", macros_printer) 
        # PETSc location if it exists.
        self.add_path("PETSC", macros_printer)
       # TRILINOS location if it exists.
        self.add_path("TRILINOS", macros_printer)
        # ALBANY location if it exists.
        self.add_path("ALBANY", macros_printer)


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
                macros_printer.print(format(model, match.text))

            # Now the same for debug mode.
            debug_matches = chain.from_iterable(
                all_matches(self.compiler_xml_tree, path, debug_dict)
                for path in paths
                )
            for match in debug_matches:
                macros_printer.print(format(model+"_DEBUG", match.text))

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
        macros_printer.print_header("External library flags.")
        add_formatted_flags(
            "SLIBS",
            lambda b, m: "add_flags(CMAKE_EXE_LINKER_FLAGS_"+b+" "+m+")"
            )
