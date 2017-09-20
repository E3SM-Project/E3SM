"""
Abstract class for restart tests

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *

logger = logging.getLogger(__name__)

class RestartTest(SystemTestsCompareTwo):

    def __init__(self, case,
                 separate_builds,
                 run_two_suffix = 'restart',
                 run_one_description = 'initial',
                 run_two_description = 'restart',
                 multisubmit = False):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds,
                                       run_two_suffix = run_two_suffix,
                                       run_one_description = run_one_description,
                                       run_two_description = run_two_description,
                                       multisubmit = multisubmit)


    def _case_one_setup(self):
        stop_n = self._case1.get_value("STOP_N")
        expect(stop_n >= 3,"STOP_N must be at least 3, STOP_N = {}".format(stop_n))

    def _case_two_setup(self):
        rest_n = self._case1.get_value("REST_N")
        stop_n = self._case1.get_value("STOP_N")
        stop_new = stop_n - rest_n
        expect(stop_new > 0, "ERROR: stop_n value {:d} too short {:d} {:d}".format(stop_new,stop_n,rest_n))
        # hist_n is set to the stop_n value of case1
        self._case.set_value("HIST_N", stop_n)
        self._case.set_value("STOP_N", stop_new)
        self._case.set_value("CONTINUE_RUN",True)
        self._case.set_value("REST_OPTION", "never")
