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
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)

class LII(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'interp',
                                       run_one_description = 'use_init_interp set to false',
                                       run_two_description = 'use_init_interp set to true')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_init_interp = .false.")

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_init_interp = .true.")

