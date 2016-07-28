"""
Implementation of the CIME LII test.

This is a CLM specific test:
Verifies that interpolation of initial conditions onto an identical
configuration gives identical results:
(1) do a run with use_init_interp false (suffix base)
(2) do a run with use_init_interp true (suffix init_interp_on)
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils import user_nl_utils

logger = logging.getLogger(__name__)

class LII(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case)

        self.set_test_suffix('interp')
        self.set_description_first_phase('use_init_interp set to false')
        self.set_description_second_phase('use_init_interp set to true')

    def _pre_build(self):
        # TODO(wjs, 2016-07-28) Rather than requiring individual tests to call
        # save_user_nl_files, should this be done always in the base class (with
        # the behavior changed so that it saves user_nl files for ALL
        # components)?
        #
        # Pro: Requires less of tests that need it
        #
        # Con: Does unnecessary work, the main downside being that it's harder
        #      to tell what's really needed for a given test
        user_nl_utils.save_user_nl_files(caseroot = self._get_caseroot(),
                                         component = "clm")

    def _common_setup(self):
        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")

    def _setup_first_phase(self):
        user_nl_utils.append_to_saved_files(caseroot = self._get_caseroot(),
                                            component = "clm",
                                            contents = "use_init_interp = .false.")

    def _setup_second_phase(self):
        user_nl_utils.append_to_saved_files(caseroot = self._get_caseroot(),
                                            component = "clm",
                                            contents = "use_init_interp = .true.")

