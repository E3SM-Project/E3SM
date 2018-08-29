"""
common utilities for buildlib
"""

from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.utils import parse_args_and_handle_standard_logging_options, setup_standard_logging_options
import sys, os, argparse
logger = logging.getLogger(__name__)

###############################################################################
def parse_input(argv):
###############################################################################

    parser = argparse.ArgumentParser()

    setup_standard_logging_options(parser)

    parser.add_argument("caseroot", default=os.getcwd(),
                        help="Case directory")

    parser.add_argument("libroot",
                        help="root for creating the library")

    parser.add_argument("bldroot",
                        help="root for building library")

    args = parse_args_and_handle_standard_logging_options(argv, parser)

    # Some compilers have trouble with long include paths, setting
    # EXEROOT to the relative path from bldroot solves the problem
    # doing it in the environment means we don't need to change all of
    # the component buildlib scripts
    with Case(args.caseroot) as case:
        os.environ["EXEROOT"] = os.path.relpath(case.get_value("EXEROOT"), args.bldroot)


    return args.caseroot, args.libroot, args.bldroot

###############################################################################
def build_cime_component_lib(case, compname, libroot, bldroot):
###############################################################################

    cimeroot  = case.get_value("CIMEROOT")
    compclass = compname[1:]
    comp_interface = case.get_value("COMP_INTERFACE")

    with open(os.path.join(bldroot,'Filepath'), 'w') as out:
        out.write(os.path.join(case.get_value('CASEROOT'), "SourceMods",
                               "src.{}\n".format(compname)) + "\n")
        if compname.startswith('d'):
            if (comp_interface == 'nuopc'):
                out.write(os.path.join(cimeroot, "src", "components", "data_comps", "dshr_nuopc") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "data_comps", compname, comp_interface) + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "data_comps", compname) + "\n")
        elif compname.startswith('x'):
            out.write(os.path.join(cimeroot, "src", "components", "xcpl_comps", "xshare") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "xcpl_comps", "xshare", comp_interface) + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "xcpl_comps", compname, comp_interface) + "\n")
        elif compname.startswith('s'):
            out.write(os.path.join(cimeroot, "src", "components", "stub_comps", compname, comp_interface) + "\n")

    # Build the component
    run_gmake(case, compclass, libroot, bldroot)

###############################################################################
def run_gmake(case, compclass, libroot, bldroot, libname="", user_cppdefs=""):
###############################################################################

    caseroot  = case.get_value("CASEROOT")
    casetools = case.get_value("CASETOOLS")
    gmake_j   = case.get_value("GMAKE_J")
    gmake     = case.get_value("GMAKE")
    mach      = case.get_value("MACH")

    complib = ""
    if libname:
        complib  = os.path.join(libroot, "lib{}.a".format(libname))
    else:
        complib  = os.path.join(libroot, "lib{}.a".format(compclass))

    makefile = os.path.join(casetools, "Makefile")
    macfile  = os.path.join(caseroot, "Macros.{}".format(mach))

    cmd = "{} complib -j {:d} MODEL={} COMPLIB={} -f {} -C {} MACFILE={} " \
        .format(gmake, gmake_j, compclass, complib, makefile, bldroot, macfile )
    if user_cppdefs:
        cmd = cmd + "USER_CPPDEFS='{}'".format(user_cppdefs )

    _, out, _ = run_cmd(cmd, combine_output=True)
    print(out.encode('utf-8'))
