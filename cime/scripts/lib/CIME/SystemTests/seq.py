"""
CIME smoke test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case_setup import case_setup, adjust_pio_layout
from CIME.check_lockedfiles import *
import shutil

logger = logging.getLogger(__name__)

class SEQ(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCommon.__init__(self, case, expected=["TEST"])

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Build two cases.
        """
        # Build the default configuration
        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
        if sharedlib_only:
            return

        # Build the model with all components with different rootpes
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")
        shutil.move("%s/%s.exe"%(exeroot,cime_model),
                    "%s/%s.exe.SEQ1"%(exeroot,cime_model))
        any_changes = False
        machpes1 = "env_mach_pes.SEQ1.xml"
        if is_locked(machpes1):
            restore(machpes1, newname="env_mach_pes.xml")
        else:
            lock_file("env_mach_pes.xml", newname=machpes1)

        comp_classes = self._case.get_values("COMP_CLASSES")
        for comp in comp_classes:
            any_changes |= self._case.get_value("ROOTPE_%s" % comp) != 0
        if any_changes:
            for comp in comp_classes:
                self._case.set_value("ROOTPE_%s"%comp, 0)
        else:
            totalpes = self._case.get_value("TOTALPES")
            newntasks = max(1, totalpes/len(comp_classes))
            rootpe = newntasks

            for comp in comp_classes:
                # here we set the cpl to have the first 2 tasks
                # and each component to have a different ROOTPE
                if comp == "CPL":
                    self._case.set_value("NTASKS_CPL", newntasks)
                else:
                    self._case.set_value("NTASKS_%s"%comp, newntasks)
                    self._case.set_value("ROOTPE_%s"%comp, rootpe)
                    rootpe += newntasks
        self._case.flush()
        case_setup(self._case, test_mode=True, reset=True)
        self.clean_build()
        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
        shutil.move("%s/%s.exe"%(exeroot,cime_model),
                    "%s/%s.exe.SEQ2"%(exeroot,cime_model))
        lock_file("env_mach_pes.xml", newname="env_mach_pes.SEQ2.xml")

    def run_phase(self):
        # Move to config_tests.xml once that's ready.
        self._case.set_value("CONTINUE_RUN", False)
        self._case.set_value("REST_OPTION", "never")
        self._case.set_value("HIST_OPTION", "$STOP_OPTION")
        self._case.set_value("HIST_N", "$STOP_N")
        self._case.flush()

        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        exeroot     = self._case.get_value("EXEROOT")
        cime_model  = self._case.get_value("MODEL")

        #
        # do an initial run test with default layout
        #
        logger.info("doing a %d %s initial test with default layout" % (stop_n, stop_option))

        shutil.copy("%s/%s.exe.SEQ1"%(exeroot,cime_model),
                    "%s/%s.exe"%(exeroot,cime_model))
        restore("env_mach_pes.SEQ1.xml", newname="env_mach_pes.xml")

        # update the pelayout settings for this run
        self._case.read_xml()
        adjust_pio_layout(self._case, self._case.get_value("PES_PER_NODE"))

        self.run_indv()

        restore("env_mach_pes.SEQ2.xml", newname="env_mach_pes.xml")

        os.remove("%s/%s.exe"%(exeroot,cime_model))
        shutil.copy("%s/%s.exe.SEQ1"%(exeroot,cime_model),
                    "%s/%s.exe"%(exeroot,cime_model))

        logger.info("doing a second %d %s test with rootpes set to zero" % (stop_n, stop_option))
        # update the pelayout settings for this run
        self._case.read_xml()
        adjust_pio_layout(self._case, self._case.get_value("PES_PER_NODE"))

        self.run_indv(suffix="seq")
        self._component_compare_test("base", "seq")
