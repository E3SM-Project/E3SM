"""
Implementation of the CIME PEM test: Tests bfb with different MPI
processor counts

This is just like running a smoke test twice - but the pe-counts
are modified the second time.
(1) Run with pes set up out of the box (suffix base)
(2) Run with half the number of tasks (suffix modpes)
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class PEM(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = 'modpes',
                                       run_one_description = 'default pe counts',
                                       run_two_description = 'halved pe counts')

    def _case_one_setup(self):
        pass

    def _case_two_setup(self):
        for comp in self._case.get_values("COMP_CLASSES"):
            ntasks = self._case.get_value("NTASKS_{}".format(comp))
            if ( ntasks > 1 ):
                self._case.set_value("NTASKS_{}".format(comp), int(ntasks/2))
        self._case.case_setup(test_mode=True, reset=True)
