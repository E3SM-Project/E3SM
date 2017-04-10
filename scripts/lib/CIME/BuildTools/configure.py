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

import shutil
from CIME.XML.standard_module_setup import *
from CIME.utils import expect
from CIME.XML.compilers import Compilers
from CIME.XML.env_mach_specific import EnvMachSpecific

logger = logging.getLogger(__name__)

def configure(machobj, output_dir, macros_format, compiler, mpilib, debug, sysos):
    """Add Macros, Depends, and env_mach_specific files to a directory.

    Arguments:
    machobj - Machines argument for this machine.
    output_dir - Directory in which to place output.
    macros_format - Container containing the string 'Makefile' to produce
                    Makefile Macros output, and/or 'CMake' for CMake output.
    compiler - String containing the compiler vendor to configure for.
    mpilib - String containing the MPI implementation to configure for.
    debug - Boolean specifying whether debugging options are enabled.
    """
    # Macros generation.
    suffixes = {'Makefile': 'make', 'CMake': 'cmake'}
    macro_maker = Compilers(machobj, compiler=compiler, mpilib=mpilib)
    for form in macros_format:
        out_file_name = os.path.join(output_dir,"Macros."+suffixes[form])
        macro_maker.write_macros_file(macros_file=out_file_name, output_format=suffixes[form])

    _copy_depends_files(machobj.get_machine_name(), machobj.machines_dir, output_dir, compiler)
    _generate_env_mach_specific(output_dir, machobj, compiler, mpilib,
                                debug, sysos)

def _copy_depends_files(machine_name, machines_dir, output_dir, compiler):
    """
    Copy any system or compiler Depends files if they do not exist in the output directory
    If there is a match for Depends.machine_name.compiler copy that and ignore the others
    """
    dfile = os.path.join(machines_dir, "Depends.%s.%s"%(machine_name,compiler))
    outputdfile = os.path.join(output_dir, "Depends.%s.%s"%(machine_name,compiler))
    if os.path.isfile(dfile):
        if not os.path.isfile(outputdfile):
            shutil.copyfile(dfile, outputdfile)
    else:
        for dep in (machine_name, compiler):
            dfile = os.path.join(machines_dir, "Depends.%s"%dep)
            outputdfile = os.path.join(output_dir, "Depends.%s"%dep)
            if os.path.isfile(dfile) and not os.path.isfile(outputdfile):
                shutil.copyfile(dfile, outputdfile)


def _generate_env_mach_specific(output_dir, machobj, compiler, mpilib, debug, sysos):
    """
    env_mach_specific generation.
    """
    ems_path = os.path.join(output_dir, "env_mach_specific.xml")
    if os.path.exists(ems_path):
        logger.warn("%s already exists, delete to replace"%ems_path)
        return
    ems_file = EnvMachSpecific(output_dir)
    ems_file.populate(machobj)
    ems_file.write()
    for shell in ('sh', 'csh'):
        ems_file.make_env_mach_specific_file(compiler, debug, mpilib, shell)
        shell_path = os.path.join(output_dir, ".env_mach_specific." + shell)
        with open(shell_path, 'a') as shell_file:
            if shell == 'sh':
                shell_file.write("\nexport COMPILER=%s\n" % compiler)
                shell_file.write("export MPILIB=%s\n" % mpilib)
                shell_file.write("export DEBUG=%s\n" % repr(debug).upper())
                shell_file.write("export OS=%s\n" % sysos)
            else:
                shell_file.write("\nsetenv COMPILER %s\n" % compiler)
                shell_file.write("setenv MPILIB %s\n" % mpilib)
                shell_file.write("setenv DEBUG %s\n" % repr(debug).upper())
                shell_file.write("setenv OS %s\n" % sysos)
