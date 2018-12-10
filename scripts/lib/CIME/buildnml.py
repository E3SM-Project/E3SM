"""
common implementation for building namelist commands

These are used by components/<model_type>/<component>/cime_config/buildnml
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, parse_args_and_handle_standard_logging_options, setup_standard_logging_options
import sys, os, argparse

logger = logging.getLogger(__name__)

###############################################################################
def parse_input(argv):
###############################################################################

    parser = argparse.ArgumentParser()

    setup_standard_logging_options(parser)

    parser.add_argument("caseroot", default=os.getcwd(),
                        help="Case directory")

    args = parse_args_and_handle_standard_logging_options(argv, parser)

    return args.caseroot

###############################################################################
#pylint: disable=unused-argument
def build_xcpl_nml(case, caseroot, compname):
###############################################################################
    compclasses = case.get_values("COMP_CLASSES")
    compclass = None
    for compclass in compclasses:
        if case.get_value("COMP_{}".format(compclass)) == compname:
            break
    expect(compclass is not None,
           "Could not identify compclass for compname {}".format(compname))
    rundir = case.get_value("RUNDIR")
    comp_interface = case.get_value("COMP_INTERFACE")

    if comp_interface != "nuopc":
        ninst  = case.get_value("NINST_{}".format(compclass.upper()))
    else:
        ninst = case.get_value("NINST")
    if not ninst:
        ninst = 1
    nx     = case.get_value("{}_NX".format(compclass.upper()))
    ny     = case.get_value("{}_NY".format(compclass.upper()))
    if compname == "xrof":
        flood_mode = case.get_value('XROF_FLOOD_MODE')

    extras = []
    dtype = 1
    npes = 0
    length = 0

    if compname == "xatm":
        if ny == 1:
            dtype = 2
        extras = [["24",
                   "ncpl  number of communications w/coupler per dat"],
                  ["0.0",
                   "simul time proxy (secs): time between cpl comms"]]
    elif compname == "xglc" or compname == "xice":
        dtype = 2
    elif compname == "xlnd":
        dtype = 11
    elif compname == "xocn":
        dtype = 4
    elif compname == "xrof":
        dtype = 11
        if flood_mode == "ACTIVE":
            extras = [[".true.", "flood flag"]]
        else:
            extras = [[".false.", "flood flag"]]

    for i in range(1, ninst + 1):
        # If only 1 file, name is 'compclass_in'
        # otherwise files are 'compclass_in0001', 'compclass_in0002', etc
        if ninst == 1:
            filename = os.path.join(rundir, "{}_in".format(compname))
        else:
            filename = os.path.join(rundir, "{}_in_{:04d}".format(compname, i))

        with open(filename, 'w') as infile:
            infile.write("{:<20d} ! i-direction global dimension\n".format(nx))
            infile.write("{:<20d} ! j-direction global dimension\n".format(ny))
            infile.write("{:<20d} ! decomp_type  1=1d-by-lat, 2=1d-by-lon, 3=2d, 4=2d evensquare, 11=segmented\n".format(dtype))
            infile.write("{:<20d} ! num of pes for i (type 3 only)\n".format(npes))
            infile.write("{:<20d} ! length of segments (type 4 only)\n".format(length))
            for extra in extras:
                #infile.write("{:-20s} ! {}\n".format(extra[0], extra[1]))
                infile.write("{:<20s} ! {}\n".format(extra[0], extra[1]))

###############################################################################
def create_namelist_infile(case, user_nl_file, namelist_infile, infile_text=""):
###############################################################################
    lines_input = []
    if os.path.isfile(user_nl_file):
        with open(user_nl_file, "r") as file_usernl:
            lines_input = file_usernl.readlines()
    else:
        logger.warning("WARNING: No file {} found in case directory".format(user_nl_file))

    lines_output = []
    lines_output.append("&comp_inparm \n")
    if infile_text:
        lines_output.append(infile_text)
        logger.debug("file_infile {} ".format(infile_text))

    for line in lines_input:
        match1 = re.search(r"^[\&\/\!]", line)
        match2 = re.search(r"\$([\w\_])+", line)
        if match1 is None and match2 is not None:
            line = case.get_resolved_value(line)
        if match1 is None:
            lines_output.append(line)

    lines_output.append("/ \n")
    with open(namelist_infile, "w") as file_infile:
        file_infile.write("\n".join(lines_output))
