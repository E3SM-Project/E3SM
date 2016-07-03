"""
CIME ERP test.  This class inherits from SystemTestsCommon

This is a pes counts hybrid (open-MP/MPI) restart bfb test from
startup.  This is just like an ERS test but the pe-counts/threading
count are modified on retart.
(1) Do an initial run with pes set up out of the box (suffix base)
(2) Do a restart test with half the number of tasks and threads (suffix rest)
"""

import shutil, filecmp
from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.case_setup import case_setup
import CIME.utils
from system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class ERP(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, case)
        
    def build(self, sharedlib_only=False, model_only=False):
        """
        Build two cases.  Case one uses defaults, case2 uses half the number of threads
        and tasks. This test will fail for components (e.g. pop) that do not reproduce exactly
        with different numbers of mpi tasks.
        """
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        machpes1 = os.path.join("LockedFiles","env_mach_pes.ERP1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")
        else:
            logging.warn("Copying env_mach_pes.xml to %s"%(machpes1))
            shutil.copy("env_mach_pes.xml", machpes1)

        self._case.set_value("SMP_BUILD","0")
        for bld in range(1,3):
            logging.warn("Starting bld %s"%bld)

            if (bld == 2):
                # halve the number of tasks and threads
                for comp in ['ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND']:
                    ntasks    = self._case.get_value("NTASKS_%s"%comp)
                    nthreads  = self._case.get_value("NTHRDS_%s"%comp)
                    rootpe    = self._case.get_value("ROOTPE_%s"%comp)
                    if ( nthreads > 1 ):
                        self._case.set_value("BUILD_THREADED", True)
                        self._case.set_value("NTHRDS_%s"%comp, nthreads/2)
                    if ( ntasks > 1 ):
                        self._case.set_value("NTASKS_%s"%comp, ntasks/2)
                        self._case.set_value("ROOTPE_%s"%comp, rootpe/2)

            self._case.flush()
            case_setup(self._case, test_mode=True, reset=True)
            self.clean_build()
            SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)
            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.ERP%s"%(exeroot,cime_model,bld))

            shutil.copy("env_mach_pes.xml", "env_mach_pes.xml.%s"%bld )
        #
        # Because mira/cetus interprets its run script differently than
        # other systems we need to copy the original env_mach_pes.xml back
        #
        shutil.copy(machpes1,"env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml", os.path.join("LockedFiles","env_mach_pes.xml"))

    def _erp_first_phase(self):

        # Reset beginning test settings
        expect(os.path.isfile("env_mach_pes.xml.1"),
               "ERROR: env_mach_pes.xml.1 does not exist, run case.build" )

        shutil.copy("env_mach_pes.xml.1", "env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml.1", "LockedFiles/env_mach_pes.xml")

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        exefile  = "%s/%s.exe"%(exeroot,cime_model)
        exefile1 = "%s/%s.exe.ERP1"%(exeroot,cime_model)
        if (os.path.isfile(exefile)):
            os.remove(exefile)
        shutil.copy(exefile1, exefile)

        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        expect(stop_n > 0, "Bad STOP_N: %d" % stop_n)

        rest_n = stop_n/2 + 1
        self._case.set_value("REST_N", rest_n)
        self._case.set_value("REST_OPTION", stop_option)
        self._case.set_value("HIST_N", stop_n)
        self._case.set_value("HIST_OPTION", stop_option)
        self._case.set_value("CONTINUE_RUN", False)
        self._case.flush()

        expect(stop_n > 2, "ERROR: stop_n value %d too short"%stop_n)
        logger.info("doing an %s %s initial test with restart file at %s %s"
                    %(str(stop_n), stop_option, str(rest_n), stop_option))

        return SystemTestsCommon.run(self)

    def _erp_second_phase(self):

        expect(os.path.isfile("env_mach_pes.xml.2"),
               "ERROR: env_mach_pes.xml.2 does not exist, run case.build" )
        
        shutil.copy("env_mach_pes.xml.2", "env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml.2", "LockedFiles/env_mach_pes.xml")

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        exefile  = "%s/%s.exe"%(exeroot,cime_model)
        exefile2 = "%s/%s.exe.ERP2"%(exeroot,cime_model)
        if (os.path.isfile(exefile)):
            os.remove(exefile)
        shutil.copy(exefile2, exefile)

        #FIXME - this is where env_mach_pes.xml.1 seems to be rewritten to env_mach_pes.xml
        #case_setup(self._case, test_mode=True, reset=True)

        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")

        rest_n = stop_n/2 + 1
        stop_new = stop_n - rest_n
        expect(stop_new > 0, "ERROR: stop_n value %d too short %d %d"%(stop_new,stop_n,rest_n))

        self._case.set_value("STOP_N", stop_new)
        self._case.set_value("CONTINUE_RUN", True)
        self._case.set_value("REST_OPTION","never")
        self._case.flush()

        logger.info("doing an %s %s restart test"
                    %(str(stop_n), stop_option))
        success = SystemTestsCommon._run(self, "rest")

        # Compare restart file
        if success:
            return self._component_compare_test("base", "rest")
        else:
            return False

    def run(self):
        success = self._erp_first_phase()

        if success:
            return self._erp_second_phase()
        else:
            return False

    def report(self):
        SystemTestsCommon.report(self)
