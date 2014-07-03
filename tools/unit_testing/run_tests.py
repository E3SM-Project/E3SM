#!/usr/bin/env python

# Python 3 compatible printing in Python 2.
from __future__ import print_function

#=================================================
# Standard library modules.
#=================================================

# Use optparse, because many machines seem to have an old Python without
# the newer argparse.
# In short: ditch this for argparse when everyone has Python 2.7 or later.
from optparse import OptionParser

import os
from shutil import rmtree
import subprocess
import sys
from xml.etree.ElementTree import ElementTree

#=================================================
# Utility functions for finding files.
#=================================================

def file_exists(file_name):
    """Query for whether or not the named file exists."""
    return os.access(file_name, os.F_OK)

def search_paths(file_name, path_list):
    """Return the first path in path_list that contains the named file."""
    dir_found = None
    for path in path_list:
        if file_exists(os.path.join(path, file_name)):
            dir_found = path
            break
    assert dir_found is not None, "Could not find the file "+ \
        str(file_name)+". "+"Searched in directories: "+":".join(path_list)
    return dir_found

#=================================================
# Local modules.
#=================================================

this_script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

python_dir = search_paths(
    "printer.py",
    [os.path.abspath("."),
     os.path.abspath("python"),
     this_script_dir,
     os.path.join(this_script_dir, "python"),
     ]
    )

sys.path.append(python_dir)

from machine_setup import MachineCompilerSettings
from printer import Printer
from xml_test_list import *

#=================================================
# Parse options and set up output Printer.
#=================================================

parser = OptionParser(
    usage="%prog [options]",
    description="""Within build_directory (--build-dir), runs cmake on test
specification directories (from --test-spec-dir or --xml-test-list), then
builds and runs the tests defined via CMake."""
    )
parser.add_option(
    "--build-dir", dest="build_dir",
    help="""Directory where tests are built and run. Will be created if it
does not exist."""
    )
parser.add_option(
    "--build-type", dest="build_type", default="CESM_DEBUG",
    help="""Value defined for CMAKE_BUILD_TYPE."""
    )
parser.add_option(
    "--cesm-cmake", dest="cesm_cmake_dir",
    help="""Location of CESM CMake modules.

Usually this option is unnecessary, because this script can autodetect this
location."""
    )
parser.add_option(
    "--cesm-root", dest="cesm_root_dir",
    help="""Location of CESM root directory.

Usually this option is unnecessary, because this script can autodetect this
location. If this option is set, it changes the locations used to
autodetect all other paths."""
    )
parser.add_option(
    "--clean", dest="clean", action="store_true",
    default=False,
    help="""Clean build directory before building. Removes CMake cache and
runs "make clean"."""
    )
parser.add_option(
    "--color", dest="color", action="store_true",
    default=sys.stdout.isatty(),
    help="""Turn on colorized output."""
    )
parser.add_option(
    "--no-color", dest="color", action="store_false",
    help="""Turn off colorized output."""
    )
parser.add_option(
    "--compiler", dest="compiler", default="gnu",
    help="""Compiler vendor for build (gnu, ibm, intel, nag, or pgi).

Only used for lookup in CESM Machines files."""
    )
parser.add_option(
    "--enable-genf90", dest="genf90", action="store_true",
    default=False,
    help="""Use genf90.pl to regenerate out-of-date sources from .F90.in
templates.

Not enabled by default because it creates in-source output, and because it
requires genf90.pl to be in the user's path."""
    )
parser.add_option(
    "--machines-dir", dest="machines_dir",
    help="""Location of CESM Machines directory.

In a CESM checkout this option is unnecessary, because this script can
autodetect this location."""
    )
parser.add_option(
    "--test-spec-dir", dest="test_spec_dir",
    help="""Location where tests are specified."""
    )
parser.add_option(
    "-T", "--ctest-args", dest="ctest_args",
    help="""Additional arguments to pass to CTest."""
    )
parser.add_option(
    "--use-mpi", dest="use_mpi", action="store_true",
    default=False,
    help="""Turn on MPI support for tests."""
    )
parser.add_option(
    "-v", "--verbose", dest="verbose", action="store_true",
    default=False,
    help="""Print verbose output."""
    )
parser.add_option(
    "--xml-test-list", dest="xml_test_list",
    help="""Path to an XML file listing directories to run tests from."""
    )

(options, args) = parser.parse_args()

output = Printer(color=options.color)

if len(args) != 0:
    parser.print_help()
    error_string = "\n".join(["Unrecognized argument(s) detected:"]+
                             args)
    output.print_error(error_string)
    raise Exception("Bad command-line argument.")

if options.test_spec_dir is None and options.xml_test_list is None:
    parser.print_help()
    output.print_error(
        "You must specify either --test-spec-dir or --xml-test-list."
        )
    raise Exception("Missing required argument.")

if options.build_dir is None:
    parser.print_help()
    output.print_error(
        "You must specify a build directory with --build-dir."
        )
    raise Exception("Missing required argument.")

#=================================================
# Find directory and file paths.
#=================================================

# Search for the CESM root directory.
# First check the option. If not specified, look to see if there's a tools
# directory two levels up (just as a sanity check).
if options.cesm_root_dir is not None:
    cesm_root_dir = os.path.abspath(options.cesm_root_dir)
else:
    cesm_root_guess = os.path.join(this_script_dir, "..", "..")
    if file_exists(os.path.join(cesm_root_guess, "tools")):
        cesm_root_dir = os.path.abspath(cesm_root_guess)
    else:
        cesm_root_dir = None

# CMake modules.
if options.cesm_cmake_dir is not None:
    cesm_cmake_dir = os.path.abspath(options.cesm_cmake_dir)
else:
    cesm_cmake_guesses = [
        os.getcwd(),
        os.path.abspath("cmake"),
        this_script_dir,
        os.path.join(this_script_dir, "cmake"),
        ]

    if cesm_root_dir is not None:
        cesm_cmake_guesses.append(
            os.path.join(cesm_root_dir, "scripts", "ccsm_utils", "CMake")
            )

    cesm_cmake_dir = search_paths("CESM_utils.cmake", cesm_cmake_guesses)

# CESM Machines directory.
if options.machines_dir is not None:
    machines_dir = os.path.abspath(options.machines_dir)
else:
    machines_guesses = []

    if cesm_root_dir is not None:
        machines_guesses.append(os.path.join(cesm_root_dir, "scripts",
                                             "ccsm_utils", "Machines"))

    machines_guesses.append(os.path.abspath("Machines"))

    # We may not need Machines depending on build type and environment, so
    # don't necessarily.
    try:
        machines_dir = search_paths("config_compilers.xml",
                                    machines_guesses)
    except AssertionError:
        if options.build_type.startswith("CESM"):
            raise Exception(
                "CESM build type selected, but could not find Machines."
                )
        else:
            machines_dir = None

if machines_dir is None:
    compiler_xml = None
else:
    compiler_xml = os.path.join(machines_dir, "config_compilers.xml")
    assert file_exists(compiler_xml), "Machines directory should be "+ \
        machines_dir+" but no config_compilers.xml is there!"

# Get test specification directories from command line options.
suite_specs = []

if options.xml_test_list is not None:
    test_xml_tree = ElementTree()
    test_xml_tree.parse(options.xml_test_list)
    known_paths = {
        "here": os.path.abspath(os.path.dirname(options.xml_test_list)),
        }
    suite_specs.extend(suites_from_xml(test_xml_tree, known_paths))

if options.test_spec_dir is not None:
    suite_specs.append(
        TestSuiteSpec("__command_line_test__",
                      ["__command_line_test__"],
                      [os.path.abspath(options.test_spec_dir)])
        )

# Create build directory if necessary.
build_dir = os.path.abspath(options.build_dir)

if not file_exists(build_dir):
    os.mkdir(build_dir)

# Switch to the build directory.
os.chdir(build_dir)

#=================================================
# Set the machine/compiler specific environment.
#=================================================

if machines_dir is not None:

    mach_settings = MachineCompilerSettings(options.compiler.lower(),
                                            compiler_xml,
                                            options.use_mpi)
    mach_settings.set_compiler_env()

#=================================================
# Functions to perform various stages of build.
#=================================================

def cmake_stage(name, test_spec_dir):
    """Run cmake in the current working directory.

    Arguments:
    name - Name for output messages.
    test_spec_dir - Test specification directory to run CMake on.
    """
    # Clear CMake cache.
    if options.clean:
        pwd_contents = os.listdir(os.getcwd())
        if "CMakeCache.txt" in pwd_contents:
            os.remove("CMakeCache.txt")
        if "CMakeFiles" in pwd_contents:
            rmtree("CMakeFiles")
        if "CESM_Macros.cmake" in pwd_contents:
            os.remove("CESM_Macros.cmake")

    if not file_exists("CMakeCache.txt"):

        output.print_header("Running cmake for "+name+".")

        cmake_command = [
            "cmake",
            test_spec_dir,
            "-DCESM_CMAKE_MODULE_DIRECTORY="+cesm_cmake_dir,
            "-DCMAKE_BUILD_TYPE="+options.build_type,
            ]

        if options.verbose:
            cmake_command.append("-Wdev")

        if options.genf90:
            cmake_command.append("-DENABLE_GENF90=ON")
            if cesm_root_dir is not None:
                genf90_dir = os.path.join(
                    cesm_root_dir, "tools", "cprnc", "genf90"
                    )
                cmake_command.append("-DCMAKE_PROGRAM_PATH="+genf90_dir)

        if not options.color:
            cmake_command.append("-DUSE_COLOR=OFF")

        macros_path = os.path.abspath("CESM_Macros.cmake")

        if machines_dir is not None:
            with open(macros_path, "w") as macros_file:
                mach_settings.write_cmake_macros(macros_file)

        subprocess.check_call(cmake_command)

def make_stage(name):
    """Run make in the current working directory.

    Arguments:
    name - Name for output messages.
    """
    output.print_header("Running make for "+name+".")

    if options.clean:
        subprocess.check_call(["make","clean"])

    make_command = ["make"]

    if options.verbose:
        make_command.append("VERBOSE=1")

    subprocess.check_call(make_command)

#=================================================
# Iterate over input suite specs, building the tests.
#=================================================

for spec in suite_specs:

    if not file_exists(spec.name):
        os.mkdir(spec.name)

    os.chdir(spec.name)

    for label, directory in spec:

        if not file_exists(label):
            os.mkdir(label)

        os.chdir(label)

        name = spec.name+"/"+label

        cmake_stage(name, directory)
        make_stage(name)

        os.chdir("..")

    os.chdir(build_dir)

#=================================================
# Run tests.
#=================================================

for spec in suite_specs:

    os.chdir(spec.name)

    for label, directory in spec:

        os.chdir(label)

        name = spec.name+"/"+label

        output.print_header("Running CTest tests for "+name+".")

        ctest_command = ["ctest"]

        if options.verbose:
            ctest_command.append("-VV")

        if options.ctest_args is not None:
            ctest_command.extend(options.ctest_args.split(" "))

        subprocess.call(ctest_command)

        os.chdir("..")

    os.chdir(build_dir)
