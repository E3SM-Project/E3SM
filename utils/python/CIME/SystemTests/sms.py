"""
CIME smoke test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class SMS(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def run(self):
        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")

        self._case.flush()
        logger.info("doing an %d %s initial test, no restarts written"%(stop_n, stop_option))
        SystemTestsCommon.run(self)

    def report(self):
        SystemTestsCommon.report(self)
