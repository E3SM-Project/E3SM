"""
Implemetation of CIME MCC test: Compares ensemble methods

This does two runs: In the first we run a three member ensemble using the
 MULTI_DRIVER capability, then we run a second single instance case and compare
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)


class MCC(SystemTestsCompareTwo):

    def __init__(self, case):
        self._comp_classes = []
        self._test_instances = 3
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = 'single_instance',
                                       run_two_description = 'single instance',
                                       run_one_description = 'multi driver')

    def _case_one_setup(self):
        # The multicoupler case will increase the number of tasks by the
        # number of requested couplers.
        self._case.set_value("MULTI_DRIVER",True)
        self._case.set_value("NINST", self._test_instances)

    def _case_two_setup(self):
        self._case.set_value("NINST", 1)
