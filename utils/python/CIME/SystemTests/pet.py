"""
Implementation of the CIME PET test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
import CIME.utils
from system_tests_common import SystemTestsCommon


class PET(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, caseroot, case)


    def build(self):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        for comp in ['ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND']:
            if (self._case.get_value("NTHRDS_%s"%comp) <= 1):
                self._case.set_value("NTHRDS_%s"%comp,"2")

        self._case.flush()
        run_cmd("case.setup -clean ")
        run_cmd("case.setup")
        run_cmd('case.clean_build')
        SystemTestsCommon.build(self)

    def run(self):
        SystemTestsCommon.run(self)

    def report(self):
        SystemTestsCommon.report(self)
