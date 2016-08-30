"""
Implementation of the CIME PET test.  This class inherits from SystemTestsCommon

This is an openmp test to determine that changing thread counts does not change answers.
(1) do an initial run where all components are threaded by default (suffix: base)
(2) do another initial run with nthrds=1 for all components (suffix: single_thread)
"""

from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class PET(SystemTestsCompareTwo):

    _COMPONENT_LIST = ('ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND')

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'single_thread',
                                       run_one_description = 'default threading',
                                       run_two_description = 'threads set to 1')

    def _case_one_setup(self):
        # first make sure that all components have threaded settings
        for comp in self._COMPONENT_LIST:
            if self._case.get_value("NTHRDS_%s"%comp) <= 1:
                self._case.set_value("NTHRDS_%s"%comp, 2)

        # Need to redo case_setup because we may have changed the number of threads
        case_setup(self._case, reset=True)

    def _case_two_setup(self):
        #Do a run with all threads set to 1
        for comp in self._COMPONENT_LIST:
            self._case.set_value("NTHRDS_%s"%comp, 1)

        # Need to redo case_setup because we may have changed the number of threads
        case_setup(self._case, reset=True)
