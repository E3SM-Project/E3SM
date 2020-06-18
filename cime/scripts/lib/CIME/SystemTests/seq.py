"""
sequencing bfb test (10 day seq,conc tests)
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class SEQ(SystemTestsCompareTwo):

    def __init__(self, case):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCompareTwo.__init__(self,
                                       case,
                                       separate_builds=True,
                                       run_two_suffix="seq",
                                       run_one_description = "base",
                                       run_two_description = "sequence")

    def _case_one_setup(self):
        pass

    def _case_two_setup(self):
        comp_classes = self._case.get_values("COMP_CLASSES")
        any_changes = False
        for comp in comp_classes:
            any_changes |= self._case.get_value("ROOTPE_{}".format(comp)) != 0
        if any_changes:
            for comp in comp_classes:
                self._case.set_value("ROOTPE_{}".format(comp), 0)
        else:
            totalpes = self._case.get_value("TOTALPES")
            newntasks = max(1, totalpes//len(comp_classes))
            rootpe = newntasks

            for comp in comp_classes:
                # here we set the cpl to have the first 2 tasks
                # and each component to have a different ROOTPE
                if comp == "CPL":
                    self._case.set_value("NTASKS_CPL", newntasks)
                else:
                    self._case.set_value("NTASKS_{}".format(comp), newntasks)
                    self._case.set_value("ROOTPE_{}".format(comp), rootpe)
                    rootpe += newntasks

        self._case.flush()
        self._case.case_setup(test_mode=True, reset=True)
