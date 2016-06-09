"""
common utilities for buildlib
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, run_cmd, append_status, handle_standard_logging_options
from CIME.case  import Case
import sys, os, filecmp, shutil

###############################################################################
def parse_command_line(args, parser):
###############################################################################

    parser.add_argument("caseroot", default=os.getcwd(),
                        help="Case directory")

    parser.add_argument("bldroot", 
                        help="root for building library")

    parser.add_argument("libroot", 
                        help="root for creating the library")

    args = parser.parse_args()

    handle_standard_logging_options(args)

    return args.caseroot, args.bldroot, args.libroot

###############################################################################
def build_data_model(compclass, compname, caseroot, bldroot, libroot):
###############################################################################

    case = Case(caseroot)

    CIMEROOT  = case.get_value("CIMEROOT")
    CASEBUILD = case.get_value("CASEBUILD")
    CASETOOLS = case.get_value("CASETOOLS")
    GMAKE_J   = case.get_value("GMAKE_J")
    GMAKE     = case.get_value("GMAKE")
    MACH      = case.get_value("MACH") 

    # Write directory list (Filepath)
    with open('Filepath', 'w') as out:
        out.write(os.path.join(caseroot, "SourceMods", "src.%s" %compname) + "\n")
        out.write(os.path.join(CIMEROOT, "components", "data_comps", compname) + "\n")

    # Build the component
    complib  = os.path.join(libroot, "lib%s.a" % compclass)
    makefile = os.path.join(CASETOOLS, "Makefile")
    macfile  = os.path.join(caseroot, "Macros.%s" % MACH)

    cmd = GMAKE + " complib  -j " + str(GMAKE_J) + " MODEL=" + compclass + " COMPLIB=" + complib 
    cmd = cmd + " -f " + makefile + " MACFILE=" + macfile
    
    run_cmd(cmd)

