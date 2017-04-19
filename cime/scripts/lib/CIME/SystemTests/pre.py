"""
Implementation of the CIME pause/resume test: Tests having driver
'pause' (write cpl restart file) and 'resume' (read cpl restart file)
without changing restart file. Compare to non-pause/resume run.

"""

import os.path
import logging
import glob

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.utils import expect
from CIME.hist_utils import cprnc
from CIME.case_setup import case_setup

###############################################################################
class PRE(SystemTestsCompareTwo):
###############################################################################
    """
    Implementation of the CIME pause/resume test: Tests having driver
    'pause' (write cpl restart file) and 'resume' (read cpl restart file)
    without changing restart file. Compare to non-pause/resume run.
    """

    ###########################################################################
    def __init__(self, case):
    ###########################################################################
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds=False,
                                       run_two_suffix='pr',
                                       run_one_description='no pause/resume',
                                       run_two_description='pause/resume')

    ###########################################################################
    def _case_one_setup(self):
    ###########################################################################
        case_setup(self._case, test_mode=True, reset=True)

    ###########################################################################
    def _case_two_setup(self):
    ###########################################################################
        # Set up a pause/resume run
        stopopt = self._case.get_value("STOP_OPTION")
        stopn = self._case.get_value("STOP_N")
        expect((stopn % 5) == 0, "ERROR: PRE test requires that STOP_N be divisible by five")
        pausen = stopn / 5
        expect(pausen > 0, "ERROR: pause_n (%d) must be > 0, stop_n is %d"%(pausen, stopn))
        self._case.set_value("PAUSE_OPTION", stopopt)
        self._case.set_value("PAUSE_N", pausen)
        pcl = self._case.get_value("PAUSE_COMPONENT_LIST")
        expect(pcl == "cpl", 
               "ERROR: PRE test expected PAUSE_COMPONENT_LIST = 'cpl', found '%s'"%pcl)

        self._case.flush()

    ###########################################################################
    def run_phase(self):
    ###########################################################################
        SystemTestsCompareTwo.run_phase(self)
        # Look for expected coupler restart files
        logger = logging.getLogger(__name__)
        self._activate_case1()
        rundir1 = self._case.get_value("RUNDIR")
        self._activate_case2()
        rundir2 = self._case.get_value("RUNDIR")
        should_match = (self._case.get_value("DESP_MODE") == "NOCHANGE")
        pause_comps = self._case.get_value("PAUSE_COMPONENT_LIST")
        expect((pause_comps != 'none'), "Pause/Resume (PRE) test has no pause components")
        if pause_comps == 'all':
            pause_comps = self._case.get_values("COMP_CLASSES")
        else:
            pause_comps = pause_comps.split(':')
        # End if
        for comp in pause_comps:
            restart_files_1 = glob.glob(os.path.join(rundir1, '*.%s.r.*'%comp))
            if len(restart_files_1) != 1:
                logger.error("Wrong number of case1 %s restart files, %d", comp, len(restart_files_1))
            # End if
            restart_files_2 = glob.glob(os.path.join(rundir2, '*.%s.r.*'%comp))
            if len(restart_files_2) != 5:
                logger.error("Wrong number of case2 %s restart files, %d", comp, len(restart_files_2))
            # End if
            # Do cprnc of restart files.
            rfile1 = restart_files_1[0]
            # rfile2 has to match rfile1 (same time string)
            parts = os.path.basename(rfile1).split(".")
            glob_str = "*.%s"%".".join(parts[len(parts)-4:])
            restart_files_2 = glob.glob(os.path.join(rundir2, glob_str))
            if len(restart_files_2) < 1:
                logger.error("Missing case2 restart file, %s", glob_str)
            # End if
            if len(restart_files_2) > 1:
                logger.error("Multiple case2 restart files, %s", glob_str)
            # End if
            rfile2 = restart_files_2[0]
            ok, out = cprnc(comp, rfile1, rfile2, self._case, rundir2)
            logger.warning("CPRNC result: %s, file = %s"%(ok, out))
            expect((should_match == ok),
                   "%s restart files%s match"%(comp, " do not" if should_match else ""))
        # End for
