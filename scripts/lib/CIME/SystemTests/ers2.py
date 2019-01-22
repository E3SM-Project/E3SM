"""
CIME restart test 2  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class ERS2(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the ERS2 system test
        """
        SystemTestsCommon.__init__(self, case)

    def _ers2_first_phase(self):
        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        rest_n      = self._case.get_value("REST_N")

        # Don't need restarts for first run
        self._case.set_value("REST_OPTION","never")

        expect(stop_n > 0, "Bad STOP_N: {:d}".format(stop_n))
        expect(stop_n > 2, "ERROR: stop_n value {:d} too short".format(stop_n))

        logger.info("doing an {0} {1} initial test with restart file at {2} {1}".format(str(stop_n), stop_option, str(rest_n)))
        self.run_indv()

    def _ers2_second_phase(self):
        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")

        rest_n = stop_n/2 + 1
        stop_new = rest_n

        self._case.set_value("REST_OPTION",stop_option)
        self._case.set_value("STOP_N", stop_new)
        self._case.flush()
        logger.info("doing first part {} {} restart test".format(str(stop_new), stop_option))
        self.run_indv(suffix="intermediate")

        stop_new = stop_n - rest_n
        self._case.set_value("STOP_N", stop_new)
        self._case.set_value("CONTINUE_RUN", True)
        self._case.set_value("REST_OPTION","never")

        logger.info("doing second part {} {} restart test".format(str(stop_new), stop_option))
        self.run_indv(suffix="rest")

        # Compare restart file
        self._component_compare_test("base", "rest")

    def run_phase(self):
        self._ers2_first_phase()
        self._ers2_second_phase()
