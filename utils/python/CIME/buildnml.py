from CIME.XML.standard_module_setup import *
from CIME.case import Case

class buildnml(object):
    def __init__(self, case, model, comp):
        self.case = case
        self.model = model
        self.comp = comp
        self.cimeroot = case.get_value("CIMEROOT")
        self.rundir = case.get_value("RUNDIR")
        self.ninst = case.get_value("NINST_%s" % self.comp)
        self.nx = case.get_value("%s_NX" % self.comp)
        self.ny = case.get_value("%s_NY" % self.comp)

    
class buildxnml(buildnml):
    def build_namelist(self):
        extras = []
        base_filename = os.path.join(self.rundir, "%s_in" % self.model)
        dtype = 1
        npes = 0
        length = 0
        ncpl = 24
        timeproxy = 0.0
        if self.model == "xatm":
            if self.ny == 1:
                dtype = 2
            extras = [["24",
                       "ncpl  number of communications w/coupler per dat"],
                      ["0.0",
                       "simul time proxy (secs): time between cpl comms"]]
        elif self.model == "xglc" or self.model == "xice":
            dtype = 2
        elif self.model == "xlnd":
            dtype = 11
        elif self.model == "xocn":
            dtype = 4
        elif self.model == "xrof":
            dtype = 11
            flood_mode = Case('XROF_FLOOD_MODE')
            if flood_mode == "ACTIVE":
                extras = [[".true.","flood flag"]]
            else:
                extras = [[".false.","flood flag"]]

        for i in range(1, self.ninst + 1):
            # If only 1 file, name is 'model_in'
            # otherwise files are 'model_in0001', 'model_in0002', etc
            if self.ninst == 1:
                infile = open(base_filename, 'w')
            else:
                infile = open("%s%4.4d" % (base_filename, i), 'w')
                
            infile.write("%-20d ! i-direction global dimension\n" % self.nx)
            infile.write("%-20d ! j-direction global dimension\n" % self.ny)
            infile.write("%-20d ! decomp_type  1=1d-by-lat, 2=1d-by-lon,"
                         " 3=2d, 4=2d evensquare, 11=segmented\n" % dtype)
            infile.write("%-20d ! num of pes for i (type 3 only)\n" % npes)
            infile.write("%-20d ! length of segments (type 4 only)\n" % length)
            for extra in extras:
                infile.write("%-20s ! %s\n" % (extra[0], extra[1]))
            infile.close()
