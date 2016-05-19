"""
CIME smoke test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from system_tests_common import SystemTestsCommon

import shutil

logger = logging.getLogger(__name__)

class SEQ(SystemTestsCommon):

    def __init__(self, caseroot, case):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def run(self):
        # Move to config_tests.xml once that's ready.
        self._case.set_value("CONTINUE_RUN", False)
        self._case.set_value("REST_OPTION", "never")
        self._case.set_value("HIST_OPTION", "$STOP_OPTION")
        self._case.set_value("HIST_N", "$STOP_N")
        self._case.flush()

        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")

        #
        # do an initial run test with default layout
        #
        logger.info("doing a %d %s initial test with default layout" % (stop_n, stop_option))
        success = SystemTestsCommon._run(self)

        if success:
            any_changes = False
            for comp in ['ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND', 'ESP']:
                any_changes |= self._case.get_value("ROOTPE_%s" % comp) != 0
                self._case.set_value("ROOTPE_%s" % comp, 0)

            expect(any_changes, "All ROOTPEs were already zero, we aren't testing anything")

            self._case.flush()
            shutil.copy("env_mach_pes.xml", os.path.join("LockedFiles", "env_mach_pes.xml"))

            logger.info("doing a second %d %s test with rootpes set to zero" % (stop_n, stop_option))
            success = SystemTestsCommon._run(self, "seq")

        if success:
            return self._component_compare_test("base", "seq")
        else:
            return False

    def report(self):
        SystemTestsCommon.report(self)
