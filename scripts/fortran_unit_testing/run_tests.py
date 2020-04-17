#!/usr/bin/env python
from __future__ import print_function
import os, sys
_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))
sys.path.append(os.path.join(_CIMEROOT, "scripts", "utils", "python"))
sys.path.append(os.path.join(_CIMEROOT, "scripts", "fortran_unit_testing", "python"))

from standard_script_setup import *
from CIME.BuildTools.configure import configure, FakeCase
from CIME.utils import run_cmd_no_fail, stringify_bool, expect
from CIME.XML.machines import Machines
from CIME.XML.compilers import Compilers
from CIME.XML.env_mach_specific import EnvMachSpecific
from xml_test_list import TestSuiteSpec, suites_from_xml
import socket
#=================================================
# Standard library modules.
#=================================================
from printer import Printer
from shutil import rmtree
# This violates CIME policy - move to CIME/XML directory
from xml.etree.ElementTree import ElementTree

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

    parser.add_argument("--build-optimized", action="store_true",
                        help="""By default, tests are built with debug flags.
                        If this option is provided, then tests are instead built
                        in optimized mode.""")

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
    parser.add_argument("--comp-interface",
                        default="mct",
                        help="""The cime driver/cpl interface to use."""
                        )
    parser.add_argument("--color", action="store_true",
                        default=sys.stdout.isatty(),
                        help="""Turn on colorized output."""
                        )
    parser.add_argument("--no-color",  action="store_false",
                      help="""Turn off colorized output."""
    )
    parser.add_argument("--compiler",
                        help="""Compiler vendor for build (supported depends on machine).
                        If not specified, use the default for this machine."""
                        )
    parser.add_argument("--enable-genf90", action="store_true",
                        default=True,
                        help="""Use genf90.pl to regenerate out-of-date sources from .F90.in
templates.

Not enabled by default because it creates in-source output, and because it
requires genf90.pl to be in the user's path."""
                        )

    parser.add_argument("--make-j", type=int, default=8,
                        help="""Number of processes to use for build."""
                        )

    parser.add_argument("--use-mpi", action="store_true",
                        help="""If specified, run unit tests with an mpi-enabled version
                        of pFUnit, via mpirun. (Default is to use a serial build without
                        mpirun.) This requires a pFUnit build with MPI support.""")

    parser.add_argument("--mpilib",
                        help="""MPI Library to use in build.
                        If not specified, use the default for this machine/compiler.
                        Must match an MPILIB option in config_compilers.xml.
                        e.g., for cheyenne, can use 'mpt'.
                        Only relevant if --use-mpi is specified."""
    )

    parser.add_argument("--mpirun-command",
                        help="""Command to use to run an MPI executable.
                        If not specified, uses the default for this machine.
                        Only relevant if --use-mpi is specified."""
    )
    parser.add_argument(
        "--test-spec-dir", default=".",
        help="""Location where tests are specified.
        Defaults to current directory."""
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
        help="""If specified, include OpenMP support for tests.
        (Default is to run without OpenMP support.) This requires a pFUnit build with
        OpenMP support."""
    )
    parser.add_argument(
        "--xml-test-list",
        help="""Path to an XML file listing directories to run tests from."""
        )

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args, parser)
    output = Printer(color=args.color)

    if args.xml_test_list is None and args.test_spec_dir is None:
        output.print_error(
            "You must specify either --test-spec-dir or --xml-test-list."
            )
        raise Exception("Missing required argument.")

    if args.make_j < 1:
        raise Exception("--make-j must be >= 1")

    return output, args.build_dir, args.build_optimized, args.clean,\
        args.cmake_args, args.compiler, args.enable_genf90, args.machine, args.machines_dir,\
        args.make_j, args.use_mpi, args.mpilib, args.mpirun_command, args.test_spec_dir, args.ctest_args,\
        args.use_openmp, args.xml_test_list, args.verbose, args.comp_interface


def cmake_stage(name, test_spec_dir, build_optimized, use_mpiserial, mpirun_command, output, pfunit_path,
                cmake_args=None, clean=False, verbose=False, enable_genf90=True, color=True):
    """Run cmake in the current working directory.

    Arguments:
    name - Name for output messages.
    test_spec_dir - Test specification directory to run CMake on.
    use_mpiserial (logical) - If True, we'll tell CMake to include mpi-serial for tests
                              that need it
    build_optimized (logical) - If True, we'll build in optimized rather than debug mode
    """
    if clean:
        if os.path.isfile("CMakeCache.txt"):
            os.remove("CMakeCache.txt")
        if os.path.isdir("CMakeFiles"):
            rmtree("CMakeFiles")

    if not os.path.isfile("CMakeCache.txt"):

        output.print_header("Running cmake for "+name+".")

        # This build_type only has limited uses, and should probably be removed,
        # but for now it's still needed
        if build_optimized:
            build_type = "CESM"
        else:
            build_type = "CESM_DEBUG"

        cmake_command = [
            "cmake",
            "-C Macros.cmake",
            test_spec_dir,
            "-DCIMEROOT="+_CIMEROOT,
            "-DCIME_CMAKE_MODULE_DIRECTORY="+os.path.abspath(os.path.join(_CIMEROOT,"src","CMake")),
            "-DCMAKE_BUILD_TYPE="+build_type,
            "-DPFUNIT_MPIRUN='"+mpirun_command+"'",
            "-DPFUNIT_PATH="+pfunit_path
            ]
        if use_mpiserial:
            cmake_command.append("-DUSE_MPI_SERIAL=ON")
        if verbose:
            cmake_command.append("-Wdev")

        if enable_genf90:
            cmake_command.append("-DENABLE_GENF90=ON")
            genf90_dir = os.path.join(
                _CIMEROOT,"src","externals","genf90"
                )
            cmake_command.append("-DCMAKE_PROGRAM_PATH="+genf90_dir)

        if not color:
            cmake_command.append("-DUSE_COLOR=OFF")

        if cmake_args is not None:
            cmake_command.extend(cmake_args.split(" "))

        run_cmd_no_fail(" ".join(cmake_command), combine_output=True)

def make_stage(name, output, make_j, clean=False, verbose=True):
    """Run make in the current working directory.

    Arguments:
    name - Name for output messages.
    make_j (int) - number of processes to use for make
    """
    output.print_header("Running make for "+name+".")

    if clean:
        run_cmd_no_fail("make clean")

    make_command = ["make","-j",str(make_j)]

    if verbose:
        make_command.append("VERBOSE=1")

    run_cmd_no_fail(" ".join(make_command), combine_output=True)

def find_pfunit(compilerobj, mpilib, use_openmp):
    """Find the pfunit installation we'll be using, and print its path

    Aborts if necessary information cannot be found.

    Args:
    - compilerobj: Object of type Compilers
    - mpilib: String giving the mpi library we're using
    - use_openmp: Boolean
    """
    attrs = {"MPILIB": mpilib,
             "compile_threaded": "TRUE" if use_openmp else "FALSE"
             }

    pfunit_path = compilerobj.get_optional_compiler_node("PFUNIT_PATH", attributes=attrs)
    expect(pfunit_path is not None,
           """PFUNIT_PATH not found for this machine and compiler, with MPILIB={} and compile_threaded={}.
You must specify PFUNIT_PATH in config_compilers.xml, with attributes MPILIB and compile_threaded.""".format(mpilib, attrs['compile_threaded']))
    logger.info("Using PFUNIT_PATH: {}".format(compilerobj.text(pfunit_path)))
    return compilerobj.text(pfunit_path)

#=================================================
# Iterate over input suite specs, building the tests.
#=================================================


def _main():
    output, build_dir, build_optimized, clean,\
        cmake_args, compiler, enable_genf90, machine, machines_dir,\
        make_j, use_mpi, mpilib, mpirun_command, test_spec_dir, ctest_args,\
        use_openmp, xml_test_list, verbose, comp_interface \
        = parse_command_line(sys.argv)

#=================================================
# Find directory and file paths.
#=================================================
    suite_specs = []
    # TODO: this violates cime policy of direct access to xml
    # should be moved to CIME/XML
    if xml_test_list is not None:
        test_xml_tree = ElementTree()
        test_xml_tree.parse(xml_test_list)
        known_paths = {
            "here": os.path.abspath(os.path.dirname(xml_test_list)),
            }
        suite_specs.extend(suites_from_xml(test_xml_tree, known_paths))
    if test_spec_dir is not None:
        suite_specs.append(
            TestSuiteSpec("__command_line_test__",
                          ["__command_line_test__"],
                          [os.path.abspath(test_spec_dir)])
            )


    if machines_dir is not None:
        machines_file = os.path.join(machines_dir, "config_machines.xml")
        machobj = Machines(infile=machines_file, machine=machine)
    else:
        machobj = Machines(machine=machine)

    # Create build directory if necessary.
    build_dir = os.path.abspath(build_dir)

    if not os.path.isdir(build_dir):
        os.mkdir(build_dir)

    # Switch to the build directory.
    os.chdir(build_dir)
    if clean:
        pwd_contents = os.listdir(os.getcwd())
        # Clear CMake cache.
        for file_ in pwd_contents:
            if file_ in ("Macros.cmake", "env_mach_specific.xml") \
                    or file_.startswith('Depends') or file_.startswith(".env_mach_specific"):
                os.remove(file_)

    #=================================================
    # Functions to perform various stages of build.
    #=================================================

    if not use_mpi:
        mpilib = "mpi-serial"
    elif mpilib is None:
        mpilib = machobj.get_default_MPIlib()
        logger.info("Using mpilib: {}".format(mpilib))

    if compiler is None:
        compiler = machobj.get_default_compiler()
        logger.info("Compiler is {}".format(compiler))

    compilerobj = Compilers(machobj, compiler=compiler, mpilib=mpilib)

    pfunit_path = find_pfunit(compilerobj, mpilib=mpilib, use_openmp=use_openmp)

    debug = not build_optimized
    os_ = machobj.get_value("OS")

    # Create the environment, and the Macros.cmake file
    #
    #
    configure(machobj, build_dir, ["CMake"], compiler, mpilib, debug,
              comp_interface, os_, unit_testing=True)
    machspecific = EnvMachSpecific(build_dir, unit_testing=True)

    fake_case = FakeCase(compiler, mpilib, debug, comp_interface)
    machspecific.load_env(fake_case)
    cmake_args = "{}-DOS={} -DCOMPILER={} -DDEBUG={} -DMPILIB={} -Dcompile_threaded={}".format(
        "" if not cmake_args else " ", os_, compiler, stringify_bool(debug), mpilib, stringify_bool(use_openmp))

    os.environ["UNIT_TEST_HOST"] = socket.gethostname()
    if "NETCDF_PATH" in os.environ and not "NETCDF" in os.environ:
        # The CMake Netcdf find utility that we use (from pio2) seems to key off
        # of the environment variable NETCDF, but not NETCDF_PATH
        logger.info("Setting NETCDF environment variable: {}".format(os.environ["NETCDF_PATH"]))
        os.environ["NETCDF"] = os.environ["NETCDF_PATH"]

    if "NETCDFROOT" in os.environ and not "NETCDF" in os.environ:
        # The CMake Netcdf find utility that we use (from pio2) seems to key off
        # of the environment variable NETCDF, but not NETCDFROOT
        logger.info("Setting NETCDF environment variable: {}".format(os.environ["NETCDFROOT"]))
        os.environ["NETCDF"] = os.environ["NETCDFROOT"]

    if not use_mpi:
        mpirun_command = ""
    elif mpirun_command is None:
        mpi_attribs = {
            "compiler" : compiler,
            "mpilib"   : mpilib,
            "threaded" : use_openmp,
            "comp_interface" : comp_interface,
            "unit_testing" : True
        }

        # We can get away with specifying case=None since we're using exe_only=True
        mpirun_command, _, _, _ = machspecific.get_mpirun(None, mpi_attribs, None, exe_only=True)
        mpirun_command = machspecific.get_resolved_value(mpirun_command)
        logger.info("mpirun command is '{}'".format(mpirun_command))

#=================================================
# Run tests.
#=================================================

    for spec in suite_specs:
        os.chdir(build_dir)
        if os.path.isdir(spec.name):
            if clean:
                rmtree(spec.name)

        if not os.path.isdir(spec.name):
            os.mkdir(spec.name)

        for label, directory in spec:
            os.chdir(os.path.join(build_dir,spec.name))
            if not os.path.isdir(label):
                os.mkdir(label)

            os.chdir(label)

            name = spec.name+"/"+label

            if not os.path.islink("Macros.cmake"):
                os.symlink(os.path.join(build_dir,"Macros.cmake"), "Macros.cmake")
            use_mpiserial = not use_mpi
            cmake_stage(name, directory, build_optimized, use_mpiserial, mpirun_command, output, pfunit_path, verbose=verbose,
                        enable_genf90=enable_genf90, cmake_args=cmake_args)
            make_stage(name, output, make_j, clean=clean, verbose=verbose)


    for spec in suite_specs:
        os.chdir(os.path.join(build_dir,spec.name))
        for label, directory in spec:

            name = spec.name+"/"+label

            output.print_header("Running CTest tests for "+name+".")

            ctest_command = ["ctest", "--output-on-failure"]

            if verbose:
                ctest_command.append("-VV")

            if ctest_args is not None:
                ctest_command.extend(ctest_args.split(" "))

            logger.info("Running '{}'".format(" ".join(ctest_command)))
            output = run_cmd_no_fail(" ".join(ctest_command), from_dir=label, combine_output=True)
            logger.info(output)

if __name__ == "__main__":
    _main()
