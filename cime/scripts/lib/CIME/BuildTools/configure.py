#!/usr/bin/env python

"""This script writes CIME build information to a directory.

The pieces of information that will be written include:

1. Machine-specific build settings (i.e. the "Macros" file).
2. File-specific build settings (i.e. "Depends" files).
3. Environment variable loads (i.e. the env_mach_specific files).

The .env_mach_specific.sh and .env_mach_specific.csh files are specific to a
given compiler, MPI library, and DEBUG setting. By default, these will be the
machine's default compiler, the machine's default MPI library, and FALSE,
respectively. These can be changed by setting the environment variables
COMPILER, MPILIB, and DEBUG, respectively.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, safe_copy
from CIME.XML.compilers import Compilers
from CIME.XML.env_mach_specific import EnvMachSpecific

logger = logging.getLogger(__name__)

def configure(machobj, output_dir, macros_format, compiler, mpilib, debug,
              comp_interface, sysos, unit_testing=False):
    """Add Macros, Depends, and env_mach_specific files to a directory.

    Arguments:
    machobj - Machines argument for this machine.
    output_dir - Directory in which to place output.
    macros_format - Container containing the string 'Makefile' to produce
                    Makefile Macros output, and/or 'CMake' for CMake output.
    compiler - String containing the compiler vendor to configure for.
    mpilib - String containing the MPI implementation to configure for.
    debug - Boolean specifying whether debugging options are enabled.
    unit_testing - Boolean specifying whether we're running unit tests (as
                   opposed to a system run)
    """
    # Macros generation.
    suffixes = {'Makefile': 'make', 'CMake': 'cmake'}
    macro_maker = Compilers(machobj, compiler=compiler, mpilib=mpilib)
    for form in macros_format:
        out_file_name = os.path.join(output_dir,"Macros."+suffixes[form])
        macro_maker.write_macros_file(macros_file=out_file_name, output_format=suffixes[form])

    _copy_depends_files(machobj.get_machine_name(), machobj.machines_dir, output_dir, compiler)
    _generate_env_mach_specific(output_dir, machobj, compiler, mpilib,
                                debug, comp_interface, sysos, unit_testing)

def _copy_depends_files(machine_name, machines_dir, output_dir, compiler):
    """
    Copy any system or compiler Depends files if they do not exist in the output directory
    If there is a match for Depends.machine_name.compiler copy that and ignore the others
    """
    dfile = os.path.join(machines_dir, "Depends.{}.{}".format(machine_name,compiler))
    outputdfile = os.path.join(output_dir, "Depends.{}.{}".format(machine_name,compiler))
    if os.path.isfile(dfile):
        if not os.path.isfile(outputdfile):
            safe_copy(dfile, outputdfile)
    else:
        for dep in (machine_name, compiler):
            dfile = os.path.join(machines_dir, "Depends.{}".format(dep))
            outputdfile = os.path.join(output_dir, "Depends.{}".format(dep))
            if os.path.isfile(dfile) and not os.path.isfile(outputdfile):
                safe_copy(dfile, outputdfile)

class FakeCase(object):

    def __init__(self, compiler, mpilib, debug, comp_interface):
        self._vals = {"COMPILER":compiler, "MPILIB":mpilib, "DEBUG":debug, "COMP_INTERFACE":comp_interface}

    def get_value(self, attrib):
        expect(attrib in self._vals, "FakeCase does not support getting value of '%s'" % attrib)
        return self._vals[attrib]

def _generate_env_mach_specific(output_dir, machobj, compiler, mpilib, debug,
                                comp_interface, sysos, unit_testing):
    """
    env_mach_specific generation.
    """
    ems_path = os.path.join(output_dir, "env_mach_specific.xml")
    if os.path.exists(ems_path):
        logger.warning("{} already exists, delete to replace".format(ems_path))
        return
    ems_file = EnvMachSpecific(output_dir, unit_testing=unit_testing)
    ems_file.populate(machobj)
    ems_file.write()
    fake_case = FakeCase(compiler, mpilib, debug, comp_interface)
    ems_file.load_env(fake_case)
    for shell in ('sh', 'csh'):
        ems_file.make_env_mach_specific_file(shell, fake_case)
        shell_path = os.path.join(output_dir, ".env_mach_specific." + shell)
        with open(shell_path, 'a') as shell_file:
            if shell == 'sh':
                shell_file.write("\nexport COMPILER={}\n".format(compiler))
                shell_file.write("export MPILIB={}\n".format(mpilib))
                shell_file.write("export DEBUG={}\n".format(repr(debug).upper()))
                shell_file.write("export OS={}\n".format(sysos))
            else:
                shell_file.write("\nsetenv COMPILER {}\n".format(compiler))
                shell_file.write("setenv MPILIB {}\n".format(mpilib))
                shell_file.write("setenv DEBUG {}\n".format(repr(debug).upper()))
                shell_file.write("setenv OS {}\n".format(sysos))
