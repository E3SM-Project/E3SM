"""
Implementation of the CIME NCR test.  This class inherits from SystemTestsCommon

Build two exectuables for this test, the first is a default build the
second halves the number of tasks and runs two instances for each component
Lay all of the components out concurrently
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.check_lockedfiles import *

logger = logging.getLogger(__name__)

class NCR(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        machpes1 = "env_mach_pes.NCR1.xml"
        if is_locked(machpes1):
            restore(machpes1, newname="env_mach_pes.xml")

        # Build two exectuables for this test, the first is a default build, the
        # second halves the number of tasks and runs two instances for each component
        # Lay all of the components out concurrently
        for bld in range(1,3):
            logging.warn("Starting bld {}".format(bld))
            machpes = "env_mach_pes.NCR{}.xml".format(bld)
            ntasks_sum = 0
            for comp in ['ATM','OCN','WAV','GLC','ICE','ROF','LND']:
                self._case.set_value("NINST_{}".format(comp), str(bld))
                ntasks      = self._case.get_value("NTASKS_{}".format(comp))
                if(bld == 1):
                    self._case.set_value("ROOTPE_{}".format(comp), 0)
                    if ( ntasks > 1 ):
                        self._case.set_value("NTASKS_{}".format(comp), ntasks/2)
                else:
                    self._case.set_value("ROOTPE_{}".format(comp), ntasks_sum)
                    ntasks_sum += ntasks*2
                    self._case.set_value("NTASKS_{}".format(comp), ntasks*2)
            self._case.flush()

            case_setup(self._case, test_mode=True, reset=True)
            self.clean_build()
            self.build_indv(sharedlib_only, model_only)
            shutil.move("{}/{}.exe".format(exeroot,cime_model),
                        "{}/{}.exe.NCR{}".format(exeroot,cime_model,bld))
            lock_file("env_build.xml", newname="env_build.NCR{}.xml".format(bld))
            lock_file("env_mach_pes.xml", newname=machpes)

        # Because mira/cetus interprets its run script differently than
        # other systems we need to copy the original env_mach_pes.xml back
        restore(machpes1, newname="env_mach_pes.xml")

    def run_phase(self):
        os.chdir(self._caseroot)

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        # Reset beginning test settings
        expect(is_locked("env_mach_pes.NCR1.xml"),
               "ERROR: LockedFiles/env_mach_pes.NCR1.xml does not exist\n"
               "   this would been produced in the build - must run case.test_build")

        restore("env_mach_pes.NCR1.xml", newname="env_mach_pes.xml")
        restore("env_build.NCR1.xml", newname="env_build.xml")
        shutil.copy("{}/{}.exe.NCR1".format(exeroot, cime_model),
                    "{}/{}.exe".format(exeroot, cime_model))


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
        logger.info("default: doing a {} {} with NINST1".format(stop_n, stop_option))
        self.run_indv()

        #======================================================================
        # do an initial run test with NINST 2
        # want to run on same pe counts per instance and same cpl pe count
        #======================================================================

        os.remove("{}/{}.exe".format(exeroot, cime_model))
        shutil.copy("{}/{}.exe.NCR2".format(exeroot, cime_model),
                    "{}/{}.exe".format(exeroot, cime_model))
        restore("env_build.NCR2.xml", newname="env_build.xml")

        logger.info("default: doing a {} {} with NINST2".format(stop_n, stop_option))
        self.run_indv(suffix="multiinst")

        # Compare
        self._component_compare_test("base", "multiinst")
