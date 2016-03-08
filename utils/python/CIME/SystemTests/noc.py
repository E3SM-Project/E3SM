"""
CIME NOC? test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from system_tests_common import SystemTestsCommon


class NOC(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize an object interface to the NOC system test
        """
        expectedrunvars = ["CONTINUE_RUN", "REST_OPTION", "HIST_OPTION", "HIST_N"]
        SystemTestsCommon.__init__(self, caseroot, case, expectedrunvars)

