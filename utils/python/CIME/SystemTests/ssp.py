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
from CIME.utils import run_cmd

logger = logging.getLogger(__name__)

class SSP(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SSP system test
        """
        SystemTestsCommon.__init__(self, case)

    def run(self):
        caseroot = self._case.get_value("CASEROOT")
        orig_case = self._case
        orig_casevar = orig_case.get_value("CASE")

        # clone the main case to create ref1
        clone_path = "%s.ref1" % caseroot
        if os.path.exists(clone_path):
            shutil.rmtree(clone_path)
        self._case.create_clone(clone_path, keepexe=True)

        # determine run lengths needed below
        stop_nf = self._case.get_value("STOP_N")
        stop_n1 = stop_nf / 2
        stop_n2 = stop_nf - stop_n1

        #-------------------------------------------------------------------
        # (1) do a spoinup run in the main case (short term archiving is on)
        #-------------------------------------------------------------------
        os.chdir(caseroot)
        self._set_active_case(orig_case)

        logger.info("startup: doing a %s %s 00000 seconds startup run" % (stop_n1, stop_nf))
        logger.info("  writing restarts at end of run") 
        logger.info("  short term archiving is on ")

        self._case.set_value("CLM_ACCELERATED_SPINUP", "on")
        self._case.set_value("STOP_N",stop_n1)
        self._case.set_value.flush()

        dout_sr = self._case.get_value("DOUT_S_ROOT")
        success = self._run(suffix=None,
                            coupler_log_path=os.path.join(dout_sr, "logs"),
                            st_archive=True)
        if not success:
            return False

        #-------------------------------------------------------------------
        # (2) do a hybrid, non-spinup run in the cloned ref case
        #-------------------------------------------------------------------
        os.chdir(clone_path)
        self._set_active_case(clone)

        refdate = run_cmd(r'ls -1dt %s/rest/*-00000* | head -1 | sed "s/-00000.*//" | sed "s/^.*rest\///"' % dout_sr)
        refsec = "00000"

        # obtain rpointer files and necessary restart files from short term archiving directory
        rest_path = os.path.join(dout_sr, "rest", "%s-%s" % (refdate, refsec))
        for item in glob.glob("%s/*%s*" % (rest_path, refdate)):
            os.symlink(item, os.path.join(rundir, os.path.basename(item)))
        for item in glob.glob("%s/*rpointer*" % rest_path):
            shutil.copy(item, rundir)

        clone.set_value("CLM_ACCELERATED_SPINUP", "off")
        clone.set_value("RUN_TYPE", hybrid)
        clone.set_value("GET_REFCASE", False)
        clone.set_value("RUN_REFCASE", "%s" % orig_casevar)
        clone.set_value("RUN_REFDATE", refdate) 
        clone.set_value("STOP_N", stop_n2)
        clone.set_value("DOUT_S", False)
        clone.flush()

    def report(self):
        SystemTestsCommon.report(self)
