"""
Implementation of the CIME PEA test.

Builds runs and compares a single processor mpi model to a model built using mpi-serial
(1) do a run with default mpi library (suffix base)
(2) do a run with mpi-serial (suffix mpi-serial)
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.XML.machines import Machines

logger = logging.getLogger(__name__)

class PEA(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = 'mpi-serial',
                                       run_one_description = 'default mpi library',
                                       run_two_description = 'mpi-serial')

    def _common_setup(self):
        for comp in self._case.get_values("COMP_CLASSES"):
            self._case.set_value("NTASKS_{}".format(comp), 1)
            self._case.set_value("NTHRDS_{}".format(comp), 1)
            self._case.set_value("ROOTPE_{}".format(comp), 0)

    def _case_one_setup(self):
        pass

    def _case_two_setup(self):
        mach_name = self._case.get_value("MACH")
        mach_obj = Machines(machine=mach_name)
        if mach_obj.is_valid_MPIlib("mpi-serial"):
            self._case.set_value("MPILIB","mpi-serial")
        else:
            logger.warning("mpi-serial is not supported on machine '{}', "
                           "so we have to fall back to default MPI and "
                           "therefore very little is being tested".format(mach_name))

        if os.path.isfile("Macros"):
            os.remove("Macros")
        self._case.case_setup(test_mode=True, reset=True)
