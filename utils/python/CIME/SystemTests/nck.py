"""
Implementation of the CIME NCK test: Tests multi-instance

This does two runs: In the first, we use one instance per component; in the
second, we use two instances per components. NTASKS are changed in each run so
that the number of tasks per instance is the same for both runs.

Lay all of the components out sequentially
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.case_setup import case_setup

logger = logging.getLogger(__name__)

class NCK(SystemTestsCompareTwo):

    def __init__(self, case):
        self._comp_classes = []
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = 'multiinst',
                                       run_one_description = 'one instance',
                                       run_two_description = 'two instances')

    def _common_setup(self):
        # We start by halving the number of tasks for both cases. This ensures
        # that we use the same number of tasks per instance in both cases: For
        # the two-instance case, we'll double this halved number, so you may
        # think that the halving was unnecessary; but it's needed in case the
        # original NTASKS was odd. (e.g., for NTASKS originally 15, we want to
        # use NTASKS = int(15/2) * 2 = 14 tasks for case two.)
        self._comp_classes = self._case.get_value("COMP_CLASSES").split(',')
        self._comp_classes.remove("DRV")
        for comp in self._comp_classes:
            ntasks = self._case.get_value("NTASKS_%s"%comp)
            if ( ntasks > 1 ):
                self._case.set_value("NTASKS_%s"%comp, int(ntasks/2))

    def _case_one_setup(self):
        for comp in self._comp_classes:
            self._case.set_value("NINST_%s"%comp, 1)

        case_setup(self._case, test_mode=True, reset=True)

    def _case_two_setup(self):
        for comp in self._comp_classes:
            self._case.set_value("NINST_%s"%comp, 2)
            ntasks = self._case.get_value("NTASKS_%s"%comp)
            self._case.set_value("NTASKS_%s"%comp, ntasks*2)

        case_setup(self._case, test_mode=True, reset=True)
