"""
CIME ERP test.  This class inherits from SystemTestsCommon

This is a pes counts hybrid (open-MP/MPI) restart bfb test from
startup.  This is just like an ERS test but the pe-counts/threading
count are modified on retart.
(1) Do an initial run with pes set up out of the box (suffix base)
(2) Do a restart test with half the number of tasks and threads (suffix rest)
"""

import shutil
from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class ERP(SystemTestsCommon):

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
        if sharedlib_only:
            return self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        # Make backup copies of the ORIGINAL env_mach_pes.xml and
        # env_build.xml in LockedFiles if they are not there. If there
        # are already copies there then simply copy them back to
        # have the starting env_mach_pes.xml and env_build.xml
        machpes1 = os.path.join("LockedFiles","env_mach_pes.ERP1.xml")
        envbuild1 = os.path.join("LockedFiles","env_build.ERP1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")
        else:
            shutil.copy("env_mach_pes.xml","env_mach_pes.ERP1.xml")

        if ( os.path.isfile(envbuild1) ):
            shutil.copy(envbuild1,"env_build.xml")

        # Build two executables, one using the original tasks and threads (ERP1) and
        # one using the modified tasks and threads (ERP2)
        # The reason we currently need two executables that CESM-CICE has a compile time decomposition
        # For cases where ERP works, changing this decomposition will not affect answers, but it will
        # affect the executable that is used
        self._case.set_value("SMP_BUILD","0")
        for bld in range(1,3):
            logging.warn("Starting bld %s"%bld)

            if (bld == 2):
                # halve the number of tasks and threads
                for comp in self._case.get_values("COMP_CLASSES"):
                    ntasks    = self._case.get_value("NTASKS_%s"%comp)
                    nthreads  = self._case.get_value("NTHRDS_%s"%comp)
                    rootpe    = self._case.get_value("ROOTPE_%s"%comp)
                    if ( nthreads > 1 ):
                        self._case.set_value("BUILD_THREADED", True)
                        self._case.set_value("NTHRDS_%s"%comp, nthreads/2)
                    if ( ntasks > 1 ):
                        self._case.set_value("NTASKS_%s"%comp, ntasks/2)
                        self._case.set_value("ROOTPE_%s"%comp, rootpe/2)

                # Note, some components, like CESM-CICE, have
                # decomposition information in env_build.xml
                # case_setup(self._case, test_mode=True, reset=True)that
                # needs to be regenerated for the above new tasks and thread counts
                case_setup(self._case, test_mode=True, reset=True)

            # Now rebuild the system, given updated information in env_build.xml

            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
            shutil.move("%s/%s.exe"%(exeroot,cime_model),
                        "%s/%s.ERP%s.exe"%(exeroot,cime_model,bld))

            # Make copies of the new env_mach_pes.xml and the new
            # env_build.xml to be used in the run phase
            shutil.copy("env_mach_pes.xml", os.path.join("LockedFiles","env_mach_pes.ERP%s.xml"%bld ))
            shutil.copy("env_build.xml", os.path.join("LockedFiles","env_build.ERP%s.xml"%bld ))
        #
        #

    def run_phase(self):
        # run will have values 1,2
        for run in range(1,3):

            expect(os.path.isfile(os.path.join("LockedFiles","env_mach_pes.ERP%d.xml"%run)),
                   "ERROR: LockedFiles/env_mach_pes.ERP%d.xml does not exist, run case.build"%run )

            # Use the second env_mach_pes.xml and env_build.xml files
            shutil.copy(os.path.join("LockedFiles","env_mach_pes.ERP%d.xml"%run), "env_mach_pes.xml")
            shutil.copy("env_mach_pes.xml", os.path.join("LockedFiles","env_mach_pes.xml"))
            shutil.copy(os.path.join("LockedFiles","env_build.ERP%d.xml")%run, "env_build.xml")
            shutil.copy("env_build.xml", os.path.join("LockedFiles","env_build.xml"))

            # update the case to use the new values
            self._case.read_xml()

            # Use the second executable that was created
            exeroot = self._case.get_value("EXEROOT")
            cime_model = CIME.utils.get_model()
            exefile  = os.path.join(exeroot,"%s.exe"%(cime_model))
            exefile2 = os.path.join(exeroot,"%s.ERP%d.exe"%(cime_model,run))
            if (os.path.isfile(exefile)):
                os.remove(exefile)
            shutil.copy(exefile2, exefile)

            case_setup(self._case, test_mode=True, reset=True)
            stop_n      = self._case.get_value("STOP_N")
            stop_option = self._case.get_value("STOP_OPTION")

            if run == 1:
                expect(stop_n > 2, "ERROR: stop_n value %d too short"%stop_n)
                rest_n = stop_n/2 + 1
                self._case.set_value("REST_N", rest_n)
                self._case.set_value("REST_OPTION", stop_option)
                self._case.set_value("HIST_N", stop_n)
                self._case.set_value("HIST_OPTION", stop_option)
                self._case.set_value("CONTINUE_RUN", False)
                suffix = "base"
            else:
                rest_n = stop_n/2 + 1
                stop_new = stop_n - rest_n
                expect(stop_new > 0, "ERROR: stop_n value %d too short %d %d"%(stop_new,stop_n,rest_n))
                self._case.set_value("STOP_N", stop_new)
                self._case.set_value("CONTINUE_RUN", True)
                self._case.set_value("REST_OPTION","never")
                suffix = "rest"
            self.run_indv(suffix=suffix)

        self._component_compare_test("base", "rest")
