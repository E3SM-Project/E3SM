#!/usr/bin/env python

# Python 3 compatible printing in Python 2.
from __future__ import print_function

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))
sys.path.append(os.path.join(_CIMEROOT, "scripts", "utils", "python"))

from standard_script_setup import *
from CIME.BuildTools.configure import configure
#=================================================
# Standard library modules.
#=================================================

from shutil import rmtree
import subprocess

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

def parse_command_line(args):
    """Command line argument parser for configure."""
    description="""Within build_directory (--build-dir), runs cmake on test
specification directories (from --test-spec-dir or --xml-test-list), then
builds and runs the tests defined via CMake."""
    parser = argparse.ArgumentParser(description=description)

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument("--build-dir", default=".",
                        help="""Directory where tests are built and run. Will be created if it does not exist."""
                        )
    parser.add_argument("--build-type", default="CIME_DEBUG",
                        help="""Value defined for CMAKE_BUILD_TYPE."""
                        )
    parser.add_argument("--cime-cmake",
                        help="""Location of CESM CMake modules.
Usually this option is unnecessary, because this script can autodetect this
location."""
    )

    parser.add_argument("--clean", action="store_true",
                        help="""Clean build directory before building. Removes CMake cache and
runs "make clean"."""
                        )
    parser.add_argument("--cmake-args",
                        help="""Additional arguments to pass to CMake."""
                        )
    parser.add_argument("--color", action="store_true",
                        default=sys.stdout.isatty(),
                        help="""Turn on colorized output."""
                        )
    parser.add_option("--no-color",  action="store_false",
                      help="""Turn off colorized output."""
    )
    parser.add_argument("--compiler",  default="gnu",
                        help="""Compiler vendor for build (supported depends on machine).
Only used for lookup in CIME Machines files."""
                        )
    parser.add_argument("--enable-genf90", action="store_true",
                        default=True,
                        help="""Use genf90.pl to regenerate out-of-date sources from .F90.in
templates.

Not enabled by default because it creates in-source output, and because it
requires genf90.pl to be in the user's path."""
                        )

    parser.add_argument("--machines-dir",
        help="""Location of CIME Machines directory.
In a CIME checkout this option is unnecessary, because this script can
autodetect this location."""
    )
    parser.add_argument("--mpilib",
    help="""MPI Library to use in build.
Required argument (until we can get this from config_machines)
Must match an MPILIB option in config_compilers.xml.
e.g., for yellowstone, can use 'mpich2'."""
    )
    parser.add_argument(
        "--mpirun-command", default="",
        help="""Command to use to run an MPI executable.

If not specified, does not use any mpirun prefix to run executables."""
        )
    parser.add_argument(
        "--test-spec-dir",
        help="""Location where tests are specified."""
        )
    parser.add_argument(
    "-T", "--ctest-args",
    help="""Additional arguments to pass to CTest."""
    )
    parser.add_argument(
        "--use-env-compiler",  action="store_true",
        default=False,
        help="""Always use environment settings to set compiler commands.

This is only necessary if using a CIME build type, if the user wants to
override the command provided by Machines."""
    )
    parser.add_argument(
        "--use-openmp",  action="store_true",
        default=False,
        help="""Turn on OPENMP support for tests."""
    )
    parser.add_argument(
        "--xml-test-list",
        help="""Path to an XML file listing directories to run tests from."""
        )
    args = parser.parse_args()
    CIME.utils.handle_standard_logging_options(args)
    output = Printer(color=args.color)

    if args.xml_test_list is None and args.test_spec_dir is None:
        output.print_error(
            "You must specify either --test-spec-dir or --xml-test-list."
            )
        raise Exception("Missing required argument.")

    return output, args.build_dir, args.build_type, args.cime_cmake, args.clean,\
        args.cmake_args, args.compiler, args.enable_genf90, args.machines_dir,\
        args.mpilib, args.mpirun_command, args.test_spec_dir, args.ctest_args,\
        args.use_env_compiler, args.use_openmp, args.xml_test_list


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
            "-DPFUNIT_MPIRUN="+options.mpirun_command,
            ]

        if options.verbose:
            cmake_command.append("-Wdev")

        if options.genf90:
            cmake_command.append("-DENABLE_GENF90=ON")
            if cesm_root_dir is not None:
                genf90_dir = os.path.join(
                    cesm_root_dir, "cime", "externals", "genf90"
                    )
                cmake_command.append("-DCMAKE_PROGRAM_PATH="+genf90_dir)

        if not options.color:
            cmake_command.append("-DUSE_COLOR=OFF")

        if options.cmake_args is not None:
            cmake_command.extend(options.cmake_args.split(" "))

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


def _main():
    output, build_dir, build_type, cime_cmake, clean,\
        cmake_args, compiler, enable_genf90, machines_dir,\
        mpilib, mpirun_command, test_spec_dir, ctest_args,\
        use_env_compiler, use_openmp, xml_test_list \
        = parse_command_line(sys.argv)

#=================================================
# Find directory and file paths.
#=================================================

# Search for the CESM root directory.
# First check the option. If not specified, look to see if there's a tools
# directory two levels up (just as a sanity check).
    if options.cesm_root_dir is not None:
        cesm_root_dir = os.path.abspath(options.cesm_root_dir)
    else:
        cesm_root_guess = os.path.join(this_script_dir, "..", "..","..")
        if file_exists(os.path.join(cesm_root_guess, "cime")):
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
                os.path.join(cesm_root_dir, "cime", "externals", "CMake")
                )

    cesm_cmake_dir = search_paths("CESM_utils.cmake", cesm_cmake_guesses)

# CESM Machines directory.
    if options.machines_dir is not None:
        machines_dir = os.path.abspath(options.machines_dir)
    else:
        machines_guesses = []

        if cesm_root_dir is not None:
            machines_guesses.append(os.path.join(cesm_root_dir, "cime", "cime_config", "cesm", "machines"))

            machines_guesses.append(os.path.abspath("machines"))

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
                                            mpilib=options.mpilib,
                                            use_env_compiler=options.use_env_compiler,
                                            use_openmp=options.use_openmp)
    mach_settings.set_compiler_env()

#=================================================
# Functions to perform various stages of build.
#=================================================

#=================================================
# Run tests.
#=================================================

for spec in suite_specs:

    os.chdir(spec.name)

    for label, directory in spec:

        os.chdir(label)

        name = spec.name+"/"+label

        output.print_header("Running CTest tests for "+name+".")

        ctest_command = ["ctest", "--output-on-failure"]

        if options.verbose:
            ctest_command.append("-VV")

        if options.ctest_args is not None:
            ctest_command.extend(options.ctest_args.split(" "))

        subprocess.call(ctest_command)

        os.chdir("..")

    os.chdir(build_dir)




if __name__ == "__main__":
    _main()
