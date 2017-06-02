"""
Implementation of the CIME SSP test.  This class inherits from SystemTestsCommon

This is a CLM specific test:
Verifies that spinup works correctly
this test is only valid for CLM compsets with CLM45
(1) do an initial spin test
    - set CLM_ACCELERATED_SPINUP to on
    - write restarts at the end of the run, turn on short term archiving
(2) do a hybrid non-spinup simulation run
    - start from the restart files generated in (1)
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
import shutil, glob, os

logger = logging.getLogger(__name__)

class SSP(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SSP system test
        """
        SystemTestsCommon.__init__(self, case)

    def run_phase(self):
        caseroot = self._case.get_value("CASEROOT")
        orig_case = self._case
        orig_casevar = self._case.get_value("CASE")

        # clone the main case to create ref1
        clone_path = "%s.ref1" % caseroot
        if os.path.exists(clone_path):
            shutil.rmtree(clone_path)
        clone = self._case.create_clone(clone_path, keepexe=True)

        # determine run lengths needed below
        stop_nf = self._case.get_value("STOP_N")
        stop_n1 = stop_nf / 2
        stop_n2 = stop_nf - stop_n1

        #-------------------------------------------------------------------
        # (1) do a spinup run in the main case in the cloned ref case
        #      (short term archiving is on)
        #-------------------------------------------------------------------
        os.chdir(clone_path)
        self._set_active_case(clone)

        logger.info("startup: doing a %s %s 00000 seconds startup run" % (stop_n1, stop_nf))
        logger.info("  writing restarts at end of run")
        logger.info("  short term archiving is on ")

        clone.set_value("CLM_ACCELERATED_SPINUP", "on")
        clone.set_value("RTM_MODE", "NULL")
        clone.set_value("STOP_N",stop_n1)
        clone.flush()

        dout_sr = clone.get_value("DOUT_S_ROOT")
        # No history files expected, set suffix=None to avoid compare error
        self.run_indv(suffix=None, st_archive=True)

        #-------------------------------------------------------------------
        # (2) do a hybrid, non-spinup run in orig_case
        #-------------------------------------------------------------------
        os.chdir(caseroot)
        self._set_active_case(orig_case)

        refdate = run_cmd_no_fail(r'ls -1dt %s/rest/*-00000* | head -1 | sed "s/-00000.*//" | sed "s/^.*rest\///"' % dout_sr)
        refsec = "00000"

        # obtain rpointer files and necessary restart files from short term archiving directory
        rundir = self._case.get_value("RUNDIR")

        rest_path = os.path.join(dout_sr, "rest", "%s-%s" % (refdate, refsec))

        for item in glob.glob("%s/*%s*" % (rest_path, refdate)):
            os.symlink(item, os.path.join(rundir, os.path.basename(item)))

        for item in glob.glob("%s/*rpointer*" % rest_path):
            shutil.copy(item, rundir)

        self._case.set_value("CLM_ACCELERATED_SPINUP", "off")
        self._case.set_value("RUN_TYPE", "hybrid")
        self._case.set_value("GET_REFCASE", False)
        self._case.set_value("RUN_REFCASE", "%s.ref1" % orig_casevar)
        self._case.set_value("RTM_MODE", "NULL")

        self._case.set_value("RUN_REFDATE", refdate)
        self._case.set_value("STOP_N", stop_n2)
        self._case.set_value("DOUT_S", False)
        self._case.flush()

        # do the restart run (short term archiving is off)
        self.run_indv()
