"""
Implementation of the CIME NCR test.  This class inherits from SystemTestsCommon

Build two exectuables for this test:
The first runs two instances for each component with the same total number of tasks,
and runs each of them concurrently
The second is a default build

NOTE: This is currently untested, and may not be working properly
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class NCR(SystemTestsCompareTwo):

    def __init__(self, case):
        """
        initialize an NCR test
        """
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = "singleinst",
                                       run_one_description = "two instances, each with the same number of tasks",
                                       run_two_description = "default build")

    def _comp_classes(self):
        # Return the components which we need to set things for
        # ESP cannot have more than one instance, so don't set anything for it
        comp_classes = self._case.get_values("COMP_CLASSES")
        if "CPL" in comp_classes:
            comp_classes.remove("CPL")
        if "ESP" in comp_classes:
            comp_classes.remove("ESP")
        return comp_classes

    def _common_setup(self):
        # Set the default number of tasks
        for comp in self._comp_classes():
            ntasks = self._case.get_value("NTASKS_{}".format(comp))
            if ntasks > 1:
                self._case.set_value("NTASKS_{}".format(comp), ntasks // 2)

    def _case_one_setup(self):
        # Set the number of instances, the ROOTPEs, and the number of tasks
        # This case should have twice the number of instances and half the number of tasks
        # All tasks should be running concurrently
        # Note that this case must be the multiinstance one
        # to correctly set the required number of nodes and avoid crashing
        ntasks_sum = 0

        for comp in self._comp_classes():
            self._case.set_value("NINST_{}".format(comp), str(2))
            self._case.set_value("ROOTPE_{}".format(comp), ntasks_sum)
            ntasks = self._case.get_value("NTASKS_{}".format(comp)) * 2
            ntasks_sum += ntasks
            self._case.set_value("NTASKS_{}".format(comp), ntasks)
        # test_mode must be False here so the case.test file is updated
        # This ensures that the correct number of nodes are used in case it's larger than in case 2

    def _case_two_setup(self):
        for comp in self._comp_classes():
            self._case.set_value("NINST_{}".format(comp), str(1))
            self._case.set_value("ROOTPE_{}".format(comp), 0)
