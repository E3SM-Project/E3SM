"""
CIME ERR test  This class inherits from ERS
ERR tests short term archiving and restart capabilities
"""
import glob, os
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.restart_tests import RestartTest
from CIME.utils import ls_sorted_by_mtime, safe_copy

logger = logging.getLogger(__name__)

class ERR(RestartTest):
    def __init__(self, case): # pylint: disable=super-init-not-called
        """
        initialize an object interface to the ERR system test
        """
        super(ERR, self).__init__(case,
                             separate_builds = False,
                             run_two_suffix = 'rest',
                             run_one_description = 'initial',
                             run_two_description = 'restart',
                             multisubmit = True)

    def _case_one_setup(self):
        super(ERR, self)._case_one_setup()
        self._case.set_value("DOUT_S", True)

    def _case_two_setup(self):
        super(ERR, self)._case_two_setup()
        self._case.set_value("DOUT_S", False)

    def _case_two_custom_prerun_action(self):
        dout_s_root = self._case1.get_value("DOUT_S_ROOT")
        rest_root = os.path.abspath(os.path.join(dout_s_root,"rest"))
        restart_list = ls_sorted_by_mtime(rest_root)
        expect(len(restart_list) >= 1, "No restart files found in {}".format(rest_root))
        self._case.restore_from_archive(rest_dir=
            os.path.join(rest_root, restart_list[0]))

    def _case_two_custom_postrun_action(self):
        # Link back to original case1 name
        # This is needed so that the necessary files are present for
        # baseline comparison and generation,
        # since some of them may have been moved to the archive directory
        for case_file in glob.iglob(os.path.join(self._case1.get_value("RUNDIR"),
                                                 "*.nc.{}".format(self._run_one_suffix))):
            orig_file = case_file[:-(1+len(self._run_one_suffix))]
            if not os.path.isfile(orig_file):
                safe_copy(case_file, orig_file)
