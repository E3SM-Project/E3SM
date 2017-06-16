"""
Implemetation of CIME MCC test: Compares ensemble methods

This does two runs: In the first we run a three member ensemble using the
original multi component single coupler method and in the second we use
the new multi coupler method.  We then compare results with the expectation that they are bfb
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.case_setup import case_setup

logger = logging.getLogger(__name__)


class MCC(SystemTestsCompareTwo):

    def __init__(self, case):
        self._comp_classes = []
        self._test_instances = 3
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'multicoupler',
                                       run_one_description = 'single instance',
                                       run_two_description = 'multi coupler')

    def _case_one_setup(self):
        # The multicoupler case will increase the number of tasks by the
        # number of requested couplers.
        self._case.set_value("NINST_CPL", self._test_instances)
        case_setup(self._case, test_mode=False, reset=True)

    def _case_two_setup(self):
        self._case.set_value("NINST_CPL", 1)
        case_setup(self._case, test_mode=True, reset=True)
