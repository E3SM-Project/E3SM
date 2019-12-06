"""
CIME ERP test.  This class inherits from RestartTest

This is a pes counts hybrid (open-MP/MPI) restart bfb test from
startup.  This is just like an ERS test but the pe-counts/threading
count are modified on restart.
(1) Do an initial run with pes set up out of the box (suffix base)
(2) Do a restart test with half the number of tasks and threads (suffix rest)
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.restart_tests import RestartTest

logger = logging.getLogger(__name__)

class ERP(RestartTest):

    def __init__(self, case):
        """
        initialize a test object
        """
        RestartTest.__init__(self, case,
                             separate_builds = True,
                             run_two_suffix = 'rest',
                             run_one_description = 'initial',
                             run_two_description = 'restart')

    def _case_two_setup(self):
        # halve the number of tasks and threads
        for comp in self._case.get_values("COMP_CLASSES"):
            ntasks    = self._case1.get_value("NTASKS_{}".format(comp))
            nthreads  = self._case1.get_value("NTHRDS_{}".format(comp))
            rootpe    = self._case1.get_value("ROOTPE_{}".format(comp))
            if ( nthreads > 1 ):
                self._case.set_value("NTHRDS_{}".format(comp), int(nthreads/2))
            if ( ntasks > 1 ):
                self._case.set_value("NTASKS_{}".format(comp), int(ntasks/2))
                self._case.set_value("ROOTPE_{}".format(comp), int(rootpe/2))

        RestartTest._case_two_setup(self)
        self._case.case_setup(test_mode=True, reset=True)
        # Note, some components, like CESM-CICE, have
        # decomposition information in env_build.xml that
        # needs to be regenerated for the above new tasks and thread counts

    def _case_one_custom_postrun_action(self):
        self.copy_case1_restarts_to_case2()
