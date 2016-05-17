"""
CIME restart test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import run_cmd, append_status

logger = logging.getLogger(__name__)

class ERS(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize an object interface to the ERS system test
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def run(self):
        stop_n = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        rest_n = stop_n/2 + 1
        self._case.set_value("REST_N",rest_n)
        self._case.flush()
        expect(stop_n > 2, "ERROR: stop_n value %d too short"%stop_n)
        logger.info("doing an %s %s initial test with restart file at %s %s"
                    %(str(stop_n), stop_option, str(rest_n), stop_option))
        SystemTestsCommon.run(self)

        expect(self.coupler_log_indicates_run_complete(),
               'ERROR: First run of ERS test Failed to complete')
        self._case.set_value("STOP_N", stop_n - rest_n)
        self._case.set_value("CONTINUE_RUN",True)
        self._case.set_value("REST_OPTION","never")
        logger.info("doing an %s %s restart test"
                    %(str(stop_n), stop_option))
        SystemTestsCommon.run(self)

        # Compare restart file
        cmd = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools",
                           "component_compare_test.sh")
        rc, out, err = run_cmd("%s -rundir %s -testcase %s -testcase_base %s -suffix1 base -suffix2 rest"
                               %(cmd, self._case.get_value('RUNDIR'), self._case.get_value('CASE'),
                                 self._case.get_value('CASEBASEID')), ok_to_fail=True)
        if rc == 0:
            append_status(out, sfile="TestStatus")
        else:
            append_status("Component_compare_test.sh failed out: %s\n\nerr: %s\n"%(out,err)
                          ,sfile="TestStatus.log")


    def report(self):
        SystemTestsCommon.report(self)
