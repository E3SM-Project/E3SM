"""
CIME smoke test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from system_tests_common import SystemTestsCommon


class ERR(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def run(self):
        self._case.set_value("CONTINUE_RUN","FALSE")
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")
        self._case.flush()
        SystemTestsCommon.run(self)

    def report(self):
        SystemTestsCommon.report(self)
