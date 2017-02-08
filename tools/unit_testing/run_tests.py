#!/usr/bin/env python
from __future__ import print_function
import os, sys
_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))
sys.path.append(os.path.join(_CIMEROOT, "scripts", "utils", "python"))
sys.path.append(os.path.join(_CIMEROOT, "tools", "unit_testing", "python"))

from standard_script_setup import *
from CIME.BuildTools.configure import configure
from CIME.utils import expect, get_cime_root, run_cmd
from CIME.XML.machines import Machines
from xml_test_list import TestSuiteSpec

#=================================================
# Standard library modules.
#=================================================
from printer import Printer
from shutil import rmtree

logger = logging.getLogger(__name__)

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

    parser.add_argument("--machine",
                        help="The machine to create build information for.")

    parser.add_argument("--machines-dir",
                        help="The machines directory to take build information "
                        "from. Overrides the CIME_MODEL environment variable, "
                        "and must be specified if that variable is not set.")

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
    parser.add_argument("--no-color",  action="store_false",
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

    return output, args.build_dir, args.build_type, args.clean,\
        args.cmake_args, args.compiler, args.enable_genf90, args.machine, args.machines_dir,\
        args.mpilib, args.mpirun_command, args.test_spec_dir, args.ctest_args,\
        args.use_env_compiler, args.use_openmp, args.xml_test_list


def cmake_stage(name, test_spec_dir, cimeroot, cmake_args=None, clean=False, verbose=False, genf90=True, color=True):
    """Run cmake in the current working directory.

    Arguments:
    name - Name for output messages.
    test_spec_dir - Test specification directory to run CMake on.
    """
    # Clear CMake cache.
    if clean:
        pwd_contents = os.listdir(os.getcwd())
        if "CMakeCache.txt" in pwd_contents:
            os.remove("CMakeCache.txt")
        if "CMakeFiles" in pwd_contents:
            rmtree("CMakeFiles")
        if "CESM_Macros.cmake" in pwd_contents:
            os.remove("CESM_Macros.cmake")

    if not os.path.isfile("CMakeCache.txt"):

        output.print_header("Running cmake for "+name+".")

        cmake_command = [
            "cmake",
            test_spec_dir,
            "-DCMAKE_MODULE_DIRECTORY="+os.path.join(cimeroot,"externals","CMake"),
            "-DCMAKE_BUILD_TYPE="+build_type,
            "-DPFUNIT_MPIRUN="+mpirun_command,
            ]

        if verbose:
            cmake_command.append("-Wdev")

        if genf90:
            cmake_command.append("-DENABLE_GENF90=ON")
            if cesm_root_dir is not None:
                genf90_dir = os.path.join(
                    cesm_root_dir, "cime", "externals", "genf90"
                    )
                cmake_command.append("-DCMAKE_PROGRAM_PATH="+genf90_dir)

        if not color:
            cmake_command.append("-DUSE_COLOR=OFF")

        if cmake_args is not None:
            cmake_command.extend(options.cmake_args.split(" "))

        macros_path = os.path.abspath("Macros.cmake")

        run_cmd(cmake_command)

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

    run_cmd(make_command)

#=================================================
# Iterate over input suite specs, building the tests.
#=================================================


def _main():
    output, build_dir, build_type, clean,\
        cmake_args, compiler, enable_genf90, machine, machines_dir,\
        mpilib, mpirun_command, test_spec_dir, ctest_args,\
        use_env_compiler, use_openmp, xml_test_list \
        = parse_command_line(sys.argv)

#=================================================
# Find directory and file paths.
#=================================================
    suite_specs = []
    if test_spec_dir is not None:
        suite_specs.append(
            TestSuiteSpec("__command_line_test__",
                          ["__command_line_test__"],
                          [os.path.abspath(test_spec_dir)])
            )




# Search for the CESM root directory.
# First check the option. If not specified, look to see if there's a tools
# directory two levels up (just as a sanity check).
    model_root_dir = os.path.abspath(os.path.join(get_cime_root(),os.pardir))
    if machine is None:
        machine="yellowstone"

    if machines_dir is not None:
        machines_file = os.path.join(machines_dir, "config_machines.xml")
        machobj = Machines(infile=machines_file, machine=machine)
    else:
        cime_model = CIME.utils.get_model()
        machobj = Machines(machine=machine)

# Create build directory if necessary.
    build_dir = os.path.abspath(build_dir)

    if not os.path.isdir(build_dir):
        os.mkdir(build_dir)

    # Switch to the build directory.
    os.chdir(build_dir)

#=================================================
# Functions to perform various stages of build.
#=================================================
    if mpilib is None:
        mpilib = machobj.get_default_MPIlib()
    if compiler is None:
        compiler = machobj.get_default_compiler()

    debug = True
    os_ = machobj.get_value("OS")

    # Create the environment, and the Macros.cmake file
    #
    #
    configure(machobj, build_dir, ["CMake"], compiler, mpilib, debug, os_)

#=================================================
# Run tests.
#=================================================

    for spec in suite_specs:
        if not os.path.isdir(spec.name):
            os.mkdir(spec.name)

        os.chdir(spec.name)

        for label, directory in spec:
            if not os.path.isdir(label):
                os.mkdir(label)

            os.chdir(label)

            name = spec.name+"/"+label

            output.print_header("Running CTest tests for "+name+".")

            ctest_command = ["ctest", "--output-on-failure"]

#            if verbose:
            ctest_command.append("-VV")

            if ctest_args is not None:
                ctest_command.extend(ctest_args.split(" "))

            run_cmd(ctest_command)

            os.chdir("..")

            os.chdir(build_dir)




if __name__ == "__main__":
    _main()
