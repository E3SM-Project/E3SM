"""
CIME production restart test  This class inherits from SystemTestsCommon
Exact restart from startup, default 2 month + 1 month
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class ERT(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the ERT system test
        """
        SystemTestsCommon.__init__(self, case)

    def _ert_first_phase(self):

        self._case.set_value("STOP_N", 2)
        self._case.set_value("STOP_OPTION", "nmonths")
        self._case.set_value("REST_N", 1)
        self._case.set_value("REST_OPTION", "nmonths")
        self._case.set_value("HIST_N", 1)
        self._case.set_value("HIST_OPTION", "nmonths")
        self._case.set_value("AVG_HIST_N", 1)
        self._case.set_value("AVG_HIST_OPTION", "nmonths")
        self._case.set_value("CONTINUE_RUN", False)
        self._case.flush()

        logger.info("doing a 2 month initial test with restart files at 1 month")
        self.run_indv()

    def _ert_second_phase(self):

        self._case.set_value("STOP_N", 1)
        self._case.set_value("CONTINUE_RUN", True)
        self._case.set_value("REST_OPTION","never")
        self._case.flush()

        logger.info("doing an 1 month restart test with no restart files")
        self.run_indv(suffix="rest")
        # Compare restart file
        self._component_compare_test("base", "rest")

    def run_phase(self):
        self._ert_first_phase()
        self._ert_second_phase()
