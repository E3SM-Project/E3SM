"""
Implementation of the CIME NCK test.  This class inherits from SystemTestsCommon

Build two exectuables for this test, the first is a default build the
second halves the number of tasks and runs two instances for each component
Lay all of the components out sequentially
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class NCK(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        '''
        build can be called once (sharedlib_only and model_only both False)
        or twice (once with each true)
        This test requires a sharedlib build for both phases
        we must handle both cases correctly
        '''
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        if not model_only:
            machpes1 = os.path.join("LockedFiles","env_mach_pes.orig.xml")
            if os.path.isfile(machpes1):
                shutil.copy(machpes1,"env_mach_pes.xml")
            else:
                shutil.copy("env_mach_pes.xml", machpes1)

        # Build two exectuables for this test, the first is a default build, the
        # second halves the number of tasks and runs two instances for each component
        # Lay all of the components out sequentially
        for bld in range(1,3):
            logging.warn("Starting bld %s"%bld)
            machpes = os.path.join("LockedFiles","env_mach_pes.NCK%s.xml"%bld)
            if model_only:
                # This file should have been created in the sharedlib_only phase
                shutil.copy(machpes,"env_mach_pes.xml")
                self._case.read_xml()
            else:
                for comp in ['ATM','OCN','WAV','GLC','ICE','ROF','LND']:
                    self._case.set_value("NINST_%s"%comp, bld)
                    ntasks      = self._case.get_value("NTASKS_%s"%comp)
                    if(bld == 1):
                        if ( ntasks > 1 ):
                            self._case.set_value("NTASKS_%s"%comp, int(ntasks/2))
                    else:
                        self._case.set_value("NTASKS_%s"%comp, ntasks*2)
                    self._case.flush()

            case_setup(self._case, test_mode=True, reset=True)
            if not sharedlib_only:
                self.clean_build()

            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
            if not model_only:
                shutil.copy("env_mach_pes.xml", machpes)
            if not sharedlib_only:
                shutil.move("%s/%s.exe"%(exeroot,cime_model),"%s/%s.exe.NCK%s"%(exeroot,cime_model,bld))
                shutil.copy("env_build.xml",os.path.join("LockedFiles","env_build.NCK%s.xml"%bld))

        # Because mira/cetus interprets its run script differently than
        # other systems we need to copy the original env_mach_pes.xml back
#        shutil.copy(machpes1,"env_mach_pes.xml")
#        shutil.copy("env_mach_pes.xml",
#                    os.path.join("LockedFiles","env_mach_pes.xml"))

    def run_phase(self):
        os.chdir(self._caseroot)

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        # Reset beginning test settings
        expect(os.path.exists("LockedFiles/env_mach_pes.NCK1.xml"),
               "ERROR: LockedFiles/env_mach_pes.NCK1.xml does not exist\n"
               "   this would been produced in the build - must run case.test_build")

        shutil.copy("LockedFiles/env_mach_pes.NCK1.xml", "env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml", "LockedFiles/env_mach_pes.xml")
        shutil.copy("%s/%s.exe.NCK1" % (exeroot, cime_model),
                    "%s/%s.exe" % (exeroot, cime_model))
        shutil.copy("LockedFiles/env_build.NCK1.xml", "env_build.xml")
        shutil.copy("env_build.xml", "LockedFiles/env_build.xml")

        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")

        self._case.set_value("HIST_N", stop_n)
        self._case.set_value("HIST_OPTION", stop_option)
        self._case.set_value("CONTINUE_RUN", False)
        self._case.set_value("REST_OPTION", "none")
        self._case.flush()

        #======================================================================
        # do an initial run test with NINST 1
        #======================================================================
        logger.info("default: doing a %s %s with NINST1" % (stop_n, stop_option))
        self.run_indv()

        #======================================================================
        # do an initial run test with NINST 2
        # want to run on same pe counts per instance and same cpl pe count
        #======================================================================

        os.remove("%s/%s.exe" % (exeroot, cime_model))
        shutil.copy("%s/%s.exe.NCK2" % (exeroot, cime_model),
                    "%s/%s.exe" % (exeroot, cime_model))
        shutil.copy("LockedFiles/env_build.NCK2.xml", "env_build.xml")
        shutil.copy("env_build.xml", "LockedFiles/env_build.xml")
        shutil.copy("LockedFiles/env_mach_pes.NCK2.xml", "env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml", "LockedFiles/env_mach_pes.xml")

        logger.info("default: doing a %s %s with NINST2" % (stop_n, stop_option))
        self.run_indv(suffix="multiinst")
        self._component_compare_test("base", "multiinst")
