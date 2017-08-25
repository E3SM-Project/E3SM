"""
Implementation of the CIME IRT. (Interim Restart Test)
This test the model's restart capability as well as the short term archiver's interim restart capability

(1) Do a Run of length N with restart at N/2 and DOUT_S_SAVE_INTERIM_RESTART set to TRUE
(2) Archive Run using ST archive tools
(3) Recover first interim restart to the case2 run directory
(4) Start case2 from restart and run to the end of case1
(5) compare results.

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
from CIME.case_st_archive import case_st_archive, restore_from_archive
from CIME.utils import ls_sorted_by_mtime

logger = logging.getLogger(__name__)

class IRT(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds=False,
                                       run_two_suffix = 'restart',
                                       run_one_description = 'initial',
                                       run_two_description = 'restart',
                                       multisubmit = False)

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

    def _case_one_custom_postrun_action(self):
        case_st_archive(self._case)

    def _case_two_custom_prerun_action(self):
        dout_s_root = self._case1.get_value("DOUT_S_ROOT")
        restart_list = ls_sorted_by_mtime(os.path.join(dout_s_root,"rest"))
        logger.info("Restart directory list is {}".format(restart_list))
        expect(len(restart_list) >=2,"Expected at least two restart directories")
        # Get the older of the two restart directories
        restore_from_archive(self._case,
                             rest_dir=os.path.abspath(
                                 os.path.join(dout_s_root, "rest", restart_list[0])))
