"""
Solution reproducibility test - common parts shared by
different test methods. The CESM/ACME model's
multi-instance capability is used to conduct an ensemble
of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.

Different solution reproducibility test methods use
different namelist settings and postprocessing.
We will decide later whether separate test types
will be created for those different methods
(pergro, time step convergence, statistical methods)
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case_setup import case_setup
#import CIME.utils
#from CIME.check_lockedfiles import *

logger = logging.getLogger(__name__)

class SLR(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SLR test
        """
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):

        # Build exectuable with multiple instances
        # (in the development phase we use 3)
        # Lay all of the components out concurrently

        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        if not model_only:
            ninst = 3
            logging.warn("Starting to build multi-instance exe")
            for comp in ['ATM','OCN','WAV','GLC','ICE','ROF','LND']:
                ntasks = self._case.get_value("NTASKS_%s"%comp)
                self._case.set_value("ROOTPE_%s"%comp, 0)
                self._case.set_value("NINST_%s"%comp,  ninst)
                self._case.set_value("NTASKS_%s"%comp, ntasks*ninst)

            self._case.set_value("ROOTPE_CPL",0)
            self._case.set_value("NTASKS_CPL",ntasks*ninst)
            self._case.flush()

            case_setup(self._case, test_mode=False, reset=True)

        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def run_phase(self):
        self._case.set_value("STOP_N",6)

        ninst = 3
        for iinst in range(1, ninst+1):
            with open('user_nl_cam_'+str(iinst).zfill(4), 'w') as nlfile:
                 nlfile.write("avgflag_pertape = 'I' \n")
                 nlfile.write("nhtfrq = 1 \n")
                 nlfile.write("mfilt  = 1000 \n")
                 nlfile.write("iradsw = 2 \n")
                 nlfile.write("iradlw = 2 \n")
                 nlfile.write("ndens  = 1 \n")
                 nlfile.write("empty_htapes = .true. \n")

        self.run_indv()
