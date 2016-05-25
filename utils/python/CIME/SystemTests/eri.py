"""
CIME ERI test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from system_tests_common import SystemTestsCommon

class ERI(SystemTestsCommon):

    def __init__(self, caseroot=None, case=None):
        """
        initialize an object interface to the ERI system test
        """
        SystemTestsCommon.__init__(self, caseroot=caseroot, case=case)
        self._testname = "ERI"

    def run(self):
        SystemTestsCommon.run(self)
