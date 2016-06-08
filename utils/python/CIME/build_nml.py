"""
common implementation for building namelist commands

These are used by components/<model_type>/<component>/cime_config/buildnml
"""
from CIME.XML.standard_module_setup import *
from CIME.case import Case

def build_x_nml(caseroot, model, comp):
    case = Case(caseroot)
    rundir = case.get_value("RUNDIR")
    ninst = case.get_value("NINST_%s" % comp)
    nx = case.get_value("%s_NX" % comp)
    ny = case.get_value("%s_NY" % comp)
    extras = []
    dtype = 1
    npes = 0
    length = 0

    if model == "xatm":
        if ny == 1:
            dtype = 2
        extras = [["24",
                   "ncpl  number of communications w/coupler per dat"],
                  ["0.0",
                   "simul time proxy (secs): time between cpl comms"]]
    elif model == "xglc" or model == "xice":
        dtype = 2
    elif model == "xlnd":
        dtype = 11
    elif model == "xocn":
        dtype = 4
    elif model == "xrof":
        dtype = 11
        flood_mode = Case('XROF_FLOOD_MODE')
        if flood_mode == "ACTIVE":
            extras = [[".true.", "flood flag"]]
        else:
            extras = [[".false.", "flood flag"]]

    for i in range(1, ninst + 1):
        # If only 1 file, name is 'model_in'
        # otherwise files are 'model_in0001', 'model_in0002', etc
        if ninst == 1:
            filename = os.path.join(rundir, "%s_in" % model)
        else:
            filename = os.path.join(rundir, "%s_in%4.4d" % (model, i))


        with open(filename, 'w') as infile:
            infile.write("%-20d ! i-direction global dimension\n" % nx)
            infile.write("%-20d ! j-direction global dimension\n" % ny)
            infile.write("%-20d ! decomp_type  1=1d-by-lat, 2=1d-by-lon,"
                         " 3=2d, 4=2d evensquare, 11=segmented\n" % dtype)
            infile.write("%-20d ! num of pes for i (type 3 only)\n" % npes)
            infile.write("%-20d ! length of segments (type 4 only)\n"
                         % length)
            for extra in extras:
                infile.write("%-20s ! %s\n" % (extra[0], extra[1]))

