"""
Implementation of the CIME pause/resume test: Tests having driver
'pause' (write cpl restart file) and 'resume' (read cpl restart file)
possibly changing the restart file. Compared to non-pause/resume run.
Test can also be run with other component combinations.
Test requires DESP component to function correctly.
"""

import os.path
import logging
import glob

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.utils import expect
from CIME.hist_utils import cprnc

###############################################################################
class PRE(SystemTestsCompareTwo):
###############################################################################
    """
    Implementation of the CIME pause/resume test: Tests having driver
    'pause' (write cpl and/or other restart file(s)) and 'resume'
    (read cpl and/or other restart file(s)) possibly changing restart
    file. Compare to non-pause/resume run.
    """

    ###########################################################################
    def __init__(self, case):
    ###########################################################################
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds=False,
                                       run_two_suffix='pr',
                                       run_one_description='no pause/resume',
                                       run_two_description='pause/resume')
        self._stopopt  = ''
        self._stopn    = 0
        self._cprnc_exe = None
    ###########################################################################
    def _case_one_setup(self):
    ###########################################################################
        pass

    ###########################################################################
    def _case_two_setup(self):
    ###########################################################################
        # Set up a pause/resume run
        stopopt = self._case1.get_value("STOP_OPTION")
        stopn = self._case1.get_value("STOP_N")
        self._case.set_value("STOP_OPTION", stopopt)
        self._case.set_value("STOP_N", stopn)
        self._case.set_value("ESP_RUN_ON_PAUSE", "TRUE")
        if stopn > 3:
            pausen = 2
        else:
            pausen = 1
        # End if

        self._case.set_value("PAUSE_OPTION", stopopt)
        self._case.set_value("PAUSE_N", pausen)
        comps = self._case.get_values("COMP_CLASSES")
        pause_active = []
        for comp in comps:
            pause_active.append(self._case.get_value("PAUSE_ACTIVE_{}".format(comp)))

        expect(any(pause_active), "No pause_active flag is set")

        self._case.flush()

    ###########################################################################
    def run_phase(self): # pylint: disable=arguments-differ
    ###########################################################################
        self._activate_case2()
        should_match = (self._case.get_value("DESP_MODE") == "NOCHANGE")
        SystemTestsCompareTwo.run_phase(self, success_change=not should_match)
        # Look for expected coupler restart files
        logger = logging.getLogger(__name__)
        self._activate_case1()
        rundir1 = self._case.get_value("RUNDIR")
        self._cprnc_exe = self._case.get_value("CCSM_CPRNC")
        self._activate_case2()
        rundir2 = self._case.get_value("RUNDIR")
        compare_ok = True
        multi_driver = self._case.get_value("MULTI_DRIVER")
        comps = self._case.get_values("COMP_CLASSES")
        for comp in comps:
            if not self._case.get_value("PAUSE_ACTIVE_{}".format(comp)):
                continue
            if comp == "CPL":
                if multi_driver:
                    ninst = self._case.get_value("NINST_MAX")
                else:
                    ninst = 1
            else:
                ninst = self._case.get_value("NINST_{}".format(comp))

            comp_name = self._case.get_value('COMP_{}'.format(comp))
            for index in range(1,ninst+1):
                if ninst == 1:
                    rname = '*.{}.r.*'.format(comp_name)
                else:
                    rname = '*.{}_{:04d}.r.*'.format(comp_name, index)

                restart_files_1 = glob.glob(os.path.join(rundir1, rname))
                expect((len(restart_files_1) > 0), "No case1 restart files for {}".format(comp))
                restart_files_2 = glob.glob(os.path.join(rundir2, rname))
                expect((len(restart_files_2) > len(restart_files_1)),
                       "No pause (restart) files found in case2 for {}".format(comp))
                # Do cprnc of restart files.
                rfile1 = restart_files_1[len(restart_files_1) - 1]
                # rfile2 has to match rfile1 (same time string)
                parts = os.path.basename(rfile1).split(".")
                glob_str = "*.{}".format(".".join(parts[len(parts)-4:]))
                restart_files_2 = glob.glob(os.path.join(rundir2, glob_str))
                expect((len(restart_files_2) == 1),
                       "Missing case2 restart file, {}", glob_str)
                rfile2 = restart_files_2[0]
                ok = cprnc(comp, rfile1, rfile2, self._case, rundir2, cprnc_exe=self._cprnc_exe)[0]
                logger.warning("CPRNC result for {}: {}".format(os.path.basename(rfile1), "PASS" if (ok == should_match) else "FAIL"))
                compare_ok = compare_ok and (should_match == ok)

        expect(compare_ok,
               "Not all restart files {}".format("matched" if should_match else "failed to match"))
