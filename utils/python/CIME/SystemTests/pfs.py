"""
CIME performance test  This class inherits from SystemTestsCommon

20 days performance test, no restart files written
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class PFS(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the PFS system test
        """
        SystemTestsCommon.__init__(self, case)

    def run_phase(self):
        self._case.set_value("STOP_OPTION", "ndays")
        self._case.set_value("STOP_N", 20)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("CONTINUE_RUN", False)
        self._case.flush()

        logger.info("doing an 20 day initial test, no restarts written")
        self.run_indv()
