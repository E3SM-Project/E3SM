"""
Implementation of the CIME NCK test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.case_setup import case_setup
import CIME.utils
from system_tests_common import SystemTestsCommon

class NCK(SystemTestsCommon):

    def __init__(self, caseroot, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def build(self, sharedlib_only=False, model_only=False):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        machpes1 = os.path.join("LockedFiles","env_mach_pes.NCK1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")

        for bld in range(1,3):
            """
            Build two exectuables for this test, the first is a default build
            the second halves the number of tasks and runs two instances
            for each component
            """
            logging.warn("Starting bld %s"%bld)
            machpes = os.path.join("LockedFiles","env_mach_pes.NCK%s.xml"%bld)
            for comp in ['ATM','OCN','WAV','GLC','ICE','ROF','LND']:
                self._case.set_value("NINST_%s"%comp, bld)
                if(bld == 2):
                    ntasks      = self._case.get_value("NTASKS_%s"%comp)
                    rootpe      = self._case.get_value("ROOTPE_%s"%comp)
                    if ( ntasks > 1 ):
                        self._case.set_value("NTASKS_%s"%comp, ntasks/2)
                        self._case.set_value("ROOTPE_%s"%comp, rootpe/2)
            self._case.flush()

            case_setup(self._caseroot, test_mode=True, reset=True)
            self.clean_build()

            SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)
            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.NCK%s"%(exeroot,cime_model,bld))
            shutil.copy("env_build.xml",os.path.join("LockedFiles","env_build.NCK%s.xml"%bld))
            shutil.copy("env_mach_pes.xml", machpes)

        #
        # Because mira/cetus interprets its run script differently than
        # other systems we need to copy the original env_mach_pes.xml
        # back
        #
        shutil.copy(machpes1,"env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml",
                    os.path.join("LockedFiles","env_mach_pes.xml"))

    def run(self):
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

        # note - if you change the env_mach_pes.xml file - should always
        # rerun the following two case.setup commands to ensure that the right
        # settings are in the run script
        # note that the following two commands will eliminate all the batch files except
        # for the test file and copy the env_mach_pes.xml to the LockedFiles directory
        case_setup(self._caseroot, clean=True, test_mode=True)
        case_setup(self._caseroot)

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
        success = SystemTestsCommon.run(self)

        #======================================================================
        # do an initial run test with NINST 2
        # want to run on same pe counts per instance and same cpl pe count
        #======================================================================

        if success:
            shutil.copy("%s/%s.exe.NCK2" % (exeroot, cime_model),
                        "%s/%s.exe" % (exeroot, cime_model))
            shutil.copy("LockedFiles/env_build.NCK2.xml", "env_build.xml")
            shutil.copy("env_build.xml", "LockedFiles/env_build.xml")

            case_setup(self._caseroot, clean=True, test_mode=True)
            case_setup(self._caseroot)

            logger.info("default: doing a %s %s with NINST2" % (stop_n, stop_option))
            success = SystemTestsCommon._run(self, "multiinst")

        # Compare
        if success:
            return self._component_compare_test("base", "multiinst")
        else:
            return False

    def report(self):
        SystemTestsCommon.report(self)
