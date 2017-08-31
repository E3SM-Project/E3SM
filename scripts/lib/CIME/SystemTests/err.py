"""
CIME ERR test  This class inherits from ERS
ERR tests short term archiving and restart capabilities
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.restart_tests import RestartTest
from CIME.case_st_archive import restore_from_archive
from CIME.utils import ls_sorted_by_mtime

logger = logging.getLogger(__name__)

class ERR(RestartTest):

    def __init__(self, case): # pylint: disable=super-init-not-called
        """
        initialize an object interface to the ERR system test
        """
        RestartTest.__init__(self, case, # pylint: disable=non-parent-init-called
                             separate_builds = False,
                             run_two_suffix = 'rest',
                             run_one_description = 'initial',
                             run_two_description = 'restart',
                             multisubmit = True)

    def _case_two_custom_prerun_action(self):
        dout_s_root = self._case1.get_value("DOUT_S_ROOT")
        rest_root = os.path.abspath(os.path.join(dout_s_root,"rest"))
        restart_list = ls_sorted_by_mtime(rest_root)
        expect(len(restart_list) >= 1, "No restart files found in {}".format(rest_root))
        restore_from_archive(self._case, rest_dir=
            os.path.join(rest_root, restart_list[0]))
