"""
Implementation of the CIME IRT. (Interim Restart Test)
This test the model's restart capability as well as the short term archiver's interim restart capability

(1) Do a Run of length N with restart at N/2 and DOUT_S_SAVE_INTERIM_RESTART set to TRUE
(2) Archive Run using ST archive tools
(3) Recover first interim restart to the case2 run directory
(4) Start case2 from restart and run to the end of case1
(5) compare results.

"""

from CIME.SystemTests.restart_tests import RestartTest
from CIME.XML.standard_module_setup import *
from CIME.utils import ls_sorted_by_mtime

logger = logging.getLogger(__name__)

class IRT(RestartTest):

    def __init__(self, case):
        RestartTest.__init__(self, case,
                             separate_builds=False,
                             run_two_suffix = 'restart',
                             run_one_description = 'initial',
                             run_two_description = 'restart',
                             multisubmit = False)
        self._skip_pnl = False

    def _case_one_custom_postrun_action(self):
        self._case.case_st_archive()
        # Since preview namelist is run before _case_two_prerun_action, we need to do this here.
        dout_s_root = self._case1.get_value("DOUT_S_ROOT")
        restart_list = ls_sorted_by_mtime(os.path.join(dout_s_root,"rest"))
        logger.info("Restart directory list is {}".format(restart_list))
        expect(len(restart_list) >=2,"Expected at least two restart directories")
        # Get the older of the two restart directories
        self._case2.restore_from_archive(rest_dir=os.path.abspath(
            os.path.join(dout_s_root, "rest", restart_list[0])))
