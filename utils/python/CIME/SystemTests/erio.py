"""
ERIO tests restart with different PIO methods

This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon

import shutil

logger = logging.getLogger(__name__)

class ERIO(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCommon.__init__(self, case, expected=["TEST"])

        self._pio_types = self._case.get_env("run").get_valid_values("PIO_TYPENAME")
        self._stop_n = self._case.get_value("STOP_N")

    def _full_run(self, pio_type):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")

        stop_option = self._case.get_value("STOP_OPTION")
        expect(self._stop_n > 0, "Bad STOP_N: %d" % self._stop_n)

        # Move to config_tests.xml once that's ready
        rest_n = self._stop_n/2 + 1
        self._case.set_value("REST_N", rest_n)
        self._case.set_value("REST_OPTION", stop_option)
        self._case.set_value("HIST_N", self._stop_n)
        self._case.set_value("HIST_OPTION", stop_option)
        self._case.set_value("CONTINUE_RUN", False)
        self._case.flush()

        expect(self._stop_n > 2, "ERROR: stop_n value %d too short"%self._stop_n)
        logger.info("doing an %s %s initial test with restart file at %s %s with pio type %s"
                    %(str(self._stop_n), stop_option, str(rest_n), stop_option, pio_type))
        self.run_indv(suffix=pio_type)

    def _restart_run(self, pio_type, other_pio_type):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")
        stop_option = self._case.get_value("STOP_OPTION")

        rest_n = self._stop_n/2 + 1
        stop_new = self._stop_n - rest_n
        expect(stop_new > 0, "ERROR: stop_n value %d too short %d %d"%(stop_new,self._stop_n,rest_n))

        self._case.set_value("STOP_N", stop_new)
        self._case.set_value("CONTINUE_RUN", True)
        self._case.set_value("REST_OPTION","never")
        self._case.flush()
        logger.info("doing an %s %s restart test with %s against %s"
                    %(str(stop_new), stop_option, pio_type, other_pio_type))

        suffix = "%s.%s" % (other_pio_type, pio_type)
        self.run_indv(suffix=suffix)

        # Compare restart file
        self._component_compare_test(other_pio_type, suffix)

    def run_phase(self):

        for idx, pio_type1 in enumerate(self._pio_types):
            if pio_type1 != "default":
                self._case.set_value("PIO_TYPENAME", pio_type1)
                self._full_run(pio_type1)
                for pio_type2 in self._pio_types[idx+1:]:
                    if pio_type2 != "default":
                        self._case.set_value("PIO_TYPENAME", pio_type2)
                        self._restart_run(pio_type2, pio_type1)
