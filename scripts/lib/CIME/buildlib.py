"""
common utilities for buildlib
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import parse_args_and_handle_standard_logging_options, setup_standard_logging_options
import sys, os, argparse, doctest

logger = logging.getLogger(__name__)

###############################################################################
def parse_input(argv):
###############################################################################

    if "--test" in argv:
        test_results = doctest.testmod(verbose=True)
        sys.exit(1 if test_results.failed > 0 else 0)

    parser = argparse.ArgumentParser()

    setup_standard_logging_options(parser)

    parser.add_argument("caseroot", default=os.getcwd(),
                        help="Case directory")

    parser.add_argument("libroot",
                        help="root for creating the library")

    parser.add_argument("bldroot",
                        help="root for building library")

    args = parse_args_and_handle_standard_logging_options(argv, parser)

    return args.caseroot, args.libroot, args.bldroot

def build_cime_component_lib(case, compname, libroot, bldroot):
    cimeroot  = case.get_value("CIMEROOT")
    compclass = compname[1:]

    with open(os.path.join(bldroot,'Filepath'), 'w') as out:
        out.write(os.path.join(case.get_value('CASEROOT'), "SourceMods",
                               "src.{}\n".format(compname)) + "\n")
        if compname.startswith('d'):
            out.write(os.path.join(cimeroot, "src", "components", "data_comps", compname, "mct") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "data_comps", compname) + "\n")
        elif compname.startswith('x'):
            out.write(os.path.join(cimeroot, "src", "components", "xcpl_comps", "xshare") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "xcpl_comps",compname, "cpl") + "\n")
        elif compname.startswith('s'):
            out.write(os.path.join(cimeroot, "src", "components", "stub_comps", "xshare") + "\n")
            out.write(os.path.join(cimeroot, "src", "components", "stub_comps",compname, "cpl") + "\n")

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

    rc, out, err = run_cmd(cmd)
    expect(rc == 0, "Command {} failed rc={:d}\nout={}\nerr={}".format(cmd, rc, out, err))

    print "Command {} completed with output {}\nerr {}".format(cmd, out, err)
