"""
CIME restart test
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class ERS(SystemTestsCompareTwo):

    def __init__(self, case):
        """
        initialize an object interface to the ERS system test
        """
        SystemTestsCompareTwo.__init__(self, case,
                                       run_two_suffix = 'rest',
                                       run_one_description = 'initial run',
                                       run_two_description = 'restart test')

        self._stop_option = self._case.get_value("STOP_OPTION")
        self._stop_n = self._case.get_value("STOP_N")

        # Need stop_n > 2 in order to have at least one time unit both before
        # and after restart
        expect(self._stop_n > 2, "ERROR: stop_n value %d too short"%self._stop_n)

        self._rest_n = self._get_rest_n()

    def _get_rest_n(self):
        """
        Return an int giving rest_n for this test.

        Assumes that self._stop_n is already set.
        """

        rest_n = self._stop_n // 2 + 1
        expect(rest_n > 0, "ERROR: stop_n value %d generates rest_n too short: %d"
               %(self._stop_n, rest_n))
        return rest_n

    def _run_common_setup(self):
        self._case.set_value("HIST_OPTION", self._stop_option)
        self._case.set_value("HIST_N", self._stop_n)

    def _run_one_setup(self):
        # We shouldn't need to set STOP_N here, but we make this explicit for
        # symmetry with _run_two_setup
        self._case.set_value("STOP_N", self._stop_n)
        self._case.set_value("CONTINUE_RUN", False)
        self._case.set_value("REST_OPTION", self._stop_option)
        self._case.set_value("REST_N", self._rest_n)
        logger.info("doing a %d %s initial test with restart file at %d %s"
                    %(self._stop_n, self._stop_option,
                      self._rest_n, self._stop_option))

    def _run_two_setup(self):
        stop_new = self._stop_n - self._rest_n
        expect(stop_new > 0, "ERROR: stop_n value %d too short %d %d"
               %(stop_new, self._stop_n, self._rest_n))

        self._case.set_value("STOP_N", stop_new)
        self._case.set_value("CONTINUE_RUN", True)
        self._case.set_value("REST_OPTION","never")
        logger.info("doing a %d %s restart test"
                    %(stop_new, self._stop_option))
