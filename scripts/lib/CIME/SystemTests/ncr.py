"""
Implementation of the CIME NCR test.  This class inherits from SystemTestsCommon

Build two exectuables for this test:
The first is a default build
The second runs two instances for each component with the same total number of tasks,
and runs each of them concurrently
"""
from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.check_lockedfiles import *

logger = logging.getLogger(__name__)

class NCR(SystemTestsCompareTwo):

    def __init__(self, case):
        """
        initialize an NCR test
        """
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = "singleinst",
                                       run_one_description = "default build",
                                       run_two_description = ("half the number of tasks, " +
                                                              "twice the number of instances"))

    def _case_one_setup(self):
        # Set the number of instances, the ROOTPEs, and the number of tasks
        # The second case should have twice the number of instances and half the number of tasks
        # All tasks should be running concurrently
        comp_classes = self._case.get_values("COMP_CLASSES")
        # Is this correct?
        comp_classes.remove("CPL")

        ntasks_sum = 0
        # ['ATM','OCN','WAV','GLC','ICE','ROF','LND']
        for comp in comp_classes:
            self._case.set_value("NINST_{}".format(comp), str(1))
            self._case.set_value("ROOTPE_{}".format(comp), ntasks_sum)
            ntasks = self._case.get_value("NTASKS_{}".format(comp))
            ntasks_sum += ntasks * 2
            self._case.set_value("NTASKS_{}".format(comp), ntasks * 2)
        case_setup(self._case, test_mode = False, reset = True)

    def _case_two_setup(self):
        # Set the number of instances, the ROOTPEs, and the number of tasks
        # The first case should have mostly default settings;
        # though we apparently halve the number of tasks if greater than 1
        comp_classes = self._case.get_values("COMP_CLASSES")
        # Is this correct?
        comp_classes.remove("CPL")

        for comp in comp_classes:
            self._case.set_value("NINST_{}".format(comp), str(1))
            self._case.set_value("ROOTPE_{}".format(comp), 0)
            ntasks = self._case.get_value("NTASKS_{}".format(comp))
            if ntasks > 1:
                self._case.set_value("NTASKS_{}".format(comp), ntasks // 2)
        case_setup(self._case, test_mode = True, reset = True)
