from CIME.XML.standard_module_setup import *
from CIME.case import Case
"""
buildlibx - given the caseroot and component name(cmp), build the stub component
"""

class buildlib(object):
    def __init__(self, case, model, comp):
        self.case = case
        self.model = model
        self.comp = comp
        self.caseroot = case.get_value("CASEROOT")
        self.casetools = case.get_value("CASETOOLS")
        self.cimeroot = case.get_value("CIMEROOT")
        self.mach = case.get_value("MACH")
        self.objroot = case.get_value("OBJROOT")
        self.libroot = case.get_value("LIBROOT")
        self.gmake_j = case.get_value("GMAKE_J")
        self.gmake = case.get_value("GMAKE")

class buildxlib(buildlib):
    def buildlib(self):
        try:
            out = open('Filepath','w')
            msg = None
        except IOError, e:
            msg = str(e)

        expect(msg is None, msg)

        # Write directory list
        out.write(os.path.join(self.caseroot, "SourceMods", "src.%s", self.comp)
                  + "\n")
        out.write(os.path.join(self.cimeroot, "components", "xcpl_comps",
                               "xshare") + "\n")
        out.write(os.path.join(self.cimeroot, "components", "xcpl_comps",
                               self.comp, "cpl") + "\n")
        out.close()

        # generate macro values
        complib = os.path.join(self.libroot, "lib%s.a" % self.model)
        makefile = os.path.join(self.casetools, "Makefile")
        macfile = os.path.join(self.caseroot, "Macros.%s" % self.mach)

        # build 
        run_cmd("%s complib -j %d MODEL=%s COMPLIB=%s -f %s MACFILE=%s"
                % (self.gmake, self.gmake_j, self.model, complib,
                   makefile, macfile))

