"""
common utilities for buildlib
"""

from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.utils import parse_args_and_handle_standard_logging_options, setup_standard_logging_options, get_model, safe_copy
from CIME.build import get_standard_makefile_args

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
def build_cime_component_lib(case, compname, libroot, bldroot, use_old=True):
###############################################################################

    cimeroot  = case.get_value("CIMEROOT")
    casebuild = case.get_value("CASEBUILD")
    compclass = compname[1:] # This very hacky
    comp_interface = case.get_value("COMP_INTERFACE")
    confdir   = os.path.join(casebuild, "{}conf".format(compname))

    if not os.path.exists(confdir):
        os.mkdir(confdir)

    with open(os.path.join(confdir, 'Filepath'), 'w') as out:
        out.write(os.path.join(case.get_value('CASEROOT'), "SourceMods",
                               "src.{}\n".format(compname)) + "\n")
        if compname.startswith('d'):
            if (comp_interface == 'nuopc'):
                out.write(os.path.join(cimeroot, "src", "components", "data_comps_"+comp_interface, "dshr_nuopc") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "data_comps_"+comp_interface, compname, "src") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "data_comps_"+comp_interface, compname) + "\n")
        elif compname.startswith('x'):
            out.write(os.path.join(cimeroot, "src", "components", "xcpl_comps_"+comp_interface, "xshare") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "xcpl_comps_"+comp_interface, compname, "src") + "\n")
        elif compname.startswith('s'):
            out.write(os.path.join(cimeroot, "src", "components", "stub_comps_"+comp_interface, compname, "src") + "\n")

    with open(os.path.join(confdir, "CCSM_cppdefs"), "w") as out:
        out.write("")

    # Build the component
    if get_model() != "e3sm" or use_old:
        safe_copy(os.path.join(confdir, "Filepath"), bldroot)
        safe_copy(os.path.join(confdir, "CCSM_cppdefs"), bldroot)
        run_gmake(case, compclass, compname, libroot, bldroot)


###############################################################################
def run_gmake(case, compclass, compname, libroot, bldroot, libname="", user_cppdefs=""):
###############################################################################
    gmake_args = get_standard_makefile_args(case)

    gmake_j   = case.get_value("GMAKE_J")
    gmake     = case.get_value("GMAKE")

    complib = ""
    if libname:
        complib  = os.path.join(libroot, "lib{}.a".format(libname))
    else:
        complib  = os.path.join(libroot, "lib{}.a".format(compclass))

    makefile = os.path.join(case.get_value("CASETOOLS"), "Makefile")

    cmd = "{gmake} complib -j {gmake_j:d} MODEL={compclass} COMP_CLASS={compclass} COMP_NAME={compname} COMPLIB={complib} {gmake_args} -f {makefile} -C {bldroot} " \
        .format(gmake=gmake, gmake_j=gmake_j, compclass=compclass, compname=compname, complib=complib, gmake_args=gmake_args, makefile=makefile, bldroot=bldroot)
    if user_cppdefs:
        cmd = cmd + "USER_CPPDEFS='{}'".format(user_cppdefs )

    stat, out, err = run_cmd(cmd, combine_output=True)
    print(out)
    if stat:
        logger.info("buildlib stat={} err={}".format(stat,err))
        os.unlink(complib)
    return stat
