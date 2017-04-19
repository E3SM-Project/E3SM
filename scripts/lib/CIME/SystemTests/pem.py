"""
CIME PEM test.  This class inherits from SystemTestsCommon

This is a  modified pe counts mpi bfb test
This is just like running a smoke test twice - but the pe-counts
count are modified the second time.
(1) Do an initial run with pes set up out of the box (suffix base)
(2) Do an initial run with half the number of tasks (suffix modpes)
"""

from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class PEM(SystemTestsCompareTwo):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCompareTwo.__init__(self, case, True)

    def _case_one_setup(self):
        pass

    def _case_two_setup(self):
        for comp in self._case.get_values("COMP_CLASSES"):
            ntasks = self._case.get_value("NTASKS_%s"%comp)
            if ( ntasks > 1 ):
                self._case.set_value("NTASKS_%s"%comp, int(ntasks/2))
        case_setup(self._case, reset=True)
