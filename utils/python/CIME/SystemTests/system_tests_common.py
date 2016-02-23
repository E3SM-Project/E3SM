"""
Base class for CIME system tests
"""
from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.utils import run_cmd
import CIME.build as build

class SystemTestsCommon(object):
    def __init__(self, caseroot=os.getcwd(),case=None):
        """
        initialize a CIME system test object
        """
        self._caseroot = caseroot
        self._case = case

    def build(self):
#        (mfile, path, desc) = imp.find_module('case.build')
        build.case_build(self._caseroot, True)

    def run(self):
        run_cmd("case.run")
        return

    def report(self):
        return

    def check_mem_leak(self):
        rundir = self._case.get_value("RUNDIR")
        cpllogfile = min(glob.iglob(os.path.join(rundir,"cpl.log*")), key=os.path.getctime)
