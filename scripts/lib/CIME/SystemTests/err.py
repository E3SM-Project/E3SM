"""
CIME ERR test  This class inherits from ERS
ERR tests short term archiving and restart capabilities
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.SystemTests.irt import IRT

import shutil, glob

logger = logging.getLogger(__name__)

class ERR(IRT):

    def __init__(self, case): # pylint: disable=super-init-not-called
        """
        initialize an object interface to the ERR system test
        """
        SystemTestsCompareTwo.__init__(self, case, # pylint: disable=non-parent-init-called
                                       separate_builds = False,
                                       run_two_suffix = 'rest',
                                       run_one_description = 'initial',
                                       run_two_description = 'restart',
                                       multisubmit = True)
    def _case_one_custom_postrun_action(self):
        pass

    def _case_two_custom_prerun_action(self):
        dout_s_root = self._case1.get_value("DOUT_S_ROOT")
        rundir = self._case.get_value("RUNDIR")
        case = self._case.get_value("CASE")
        # First remove restart, history and rpointer files from the run directory
        for item in glob.iglob(os.path.join(rundir, "{}.*".format(case))):
            if not item.endswith("base"):
                os.remove(item)

        for item in glob.iglob(os.path.join(rundir, "rpointer.*")):
            os.remove(item)

        # Then replace them from the restart directory
        for item in glob.iglob(os.path.join(dout_s_root,"rest","*","*")):
            shutil.copy(item, rundir)
