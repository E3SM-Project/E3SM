"""
CIME PEM test.  This class inherits from SystemTestsCommon

This is a  modified pe counts mpi bfb test
This is just like running a smoke test twice - but the pe-counts
count are modified the second time.
(1) Do an initial run with pes set up out of the box (suffix base)
(2) Do an initial run with half the number of tasks (suffix modpes)
"""

import shutil
from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class PEM(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Build two cases.  Case one uses defaults, case2 uses half the number of threads
        and tasks. This test will fail for components (e.g. pop) that do not reproduce exactly
        with different numbers of mpi tasks.
        """
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        machpes1 = os.path.join("LockedFiles","env_mach_pes.PEM1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")
        else:
            logging.warn("Copying env_mach_pes.xml to %s"%(machpes1))
            shutil.copy("env_mach_pes.xml", machpes1)

        self._case.set_value("SMP_BUILD","0")
        for bld in range(1,3):
            logging.warn("Starting bld %s"%bld)

            if (bld == 2):
                # halve the number of tasks
                for comp in ['ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND']:
                    ntasks    = self._case.get_value("NTASKS_%s"%comp)
                    if ( ntasks > 1 ):
                        self._case.set_value("NTASKS_%s"%comp, ntasks/2)

            self._case.flush()
            case_setup(self._case, test_mode=True, reset=True)
            self.clean_build()
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.PEM%s"%(exeroot,cime_model,bld))

            shutil.copy("env_mach_pes.xml", "env_mach_pes.xml.%s"%bld )
        #
        # Because mira/cetus interprets its run script differently than
        # other systems we need to copy the original env_mach_pes.xml# back
        #
        shutil.copy(machpes1,"env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml", os.path.join("LockedFiles","env_mach_pes.xml"))

    def _pem_first_phase(self):

        # Reset beginning test settings
        expect(os.path.isfile("env_mach_pes.xml.1"),
               "ERROR: env_mach_pes.xml.1 does not exist, run case.build" )

        shutil.copy("env_mach_pes.xml.1", "env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml.1", "LockedFiles/env_mach_pes.xml")

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        exefile  = "%s/%s.exe"%(exeroot,cime_model)
        exefile1 = "%s/%s.exe.PEM1"%(exeroot,cime_model)
        if (os.path.isfile(exefile)):
            os.remove(exefile)
        shutil.copy(exefile1, exefile)

        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")
        stop_n = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")

        self._case.flush()
        logger.info("doing an %d %s initial test, no restarts written" % (stop_n, stop_option))

        self.run_indv()

    def _pem_second_phase(self):

        expect(os.path.isfile("env_mach_pes.xml.2"),
               "ERROR: env_mach_pes.xml.2 does not exist, run case.build" )

        shutil.copy("env_mach_pes.xml.2", "env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml.2", "LockedFiles/env_mach_pes.xml")

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        exefile  = "%s/%s.exe"%(exeroot,cime_model)
        exefile2 = "%s/%s.exe.PEM2"%(exeroot,cime_model)
        if (os.path.isfile(exefile)):
            os.remove(exefile)
        shutil.copy(exefile2, exefile)

        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")
        stop_n = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")

        self._case.flush()
        logger.info("doing an %d %s initial test, no restarts written" % (stop_n, stop_option))

        self.run_indv(suffix="modpes")
        self._component_compare_test("base", "modpes")

    def run_phase(self):
        self._pem_first_phase()
        self._pem_second_phase()
