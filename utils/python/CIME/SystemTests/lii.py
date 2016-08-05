"""
Implementation of the CIME LII test.

This is a CLM specific test:
Verifies that interpolation of initial conditions onto an identical
configuration gives identical results:
(1) do a run with use_init_interp false (suffix base)
(2) do a run with use_init_interp true (suffix init_interp_on)
"""

from CIME.SystemTests.system_tests_compare_two_clone import SystemTestsCompareTwoClone
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import save_user_nl_files, append_to_saved_files

logger = logging.getLogger(__name__)

class LII(SystemTestsCompareTwoClone):

    def __init__(self, case):
        SystemTestsCompareTwoClone.__init__(self, case,
                                            run_two_suffix = 'interp',
                                            run_one_description = 'use_init_interp set to false',
                                            run_two_description = 'use_init_interp set to true')

    def _common_setup(self):
        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")

    def _case_one_setup(self):
        # FIXME(wjs, 2016-08-05) Need to implement this method
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_init_interp = .false.")

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_init_interp = .true.")

