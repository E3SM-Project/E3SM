"""
CIME FUNIT test. This class inherits from SystemTestsCommon. It runs
the fortran unit tests; grid and compset are ignored.
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.build import post_build
from CIME.utils import append_testlog, get_cime_root
from CIME.test_status import *

logger = logging.getLogger(__name__)

class FUNIT(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the FUNIT system test
        """
        SystemTestsCommon.__init__(self, case)
        case.load_env()

    def build_phase(self, sharedlib_only=False, model_only=False):
        if not sharedlib_only:
            exeroot  = self._case.get_value("EXEROOT")
            logfile = os.path.join(exeroot, "funit.bldlog")
            with open(logfile, "w") as fd:
                fd.write("No-op\n")

            post_build(self._case, [logfile], build_complete=True)

    def get_test_spec_dir(self):
        """
        Override this to change what gets tested.
        """
        return get_cime_root()

    def run_phase(self):

        rundir   = self._case.get_value("RUNDIR")
        exeroot  = self._case.get_value("EXEROOT")
        mach     = self._case.get_value("MACH")

        log = os.path.join(rundir, "funit.log")
        if os.path.exists(log):
            os.remove(log)

        test_spec_dir = self.get_test_spec_dir()
        unit_test_tool = os.path.abspath(os.path.join(get_cime_root(),"scripts","fortran_unit_testing","run_tests.py"))
        args = "--build-dir {} --test-spec-dir {} --machine {}".format(exeroot, test_spec_dir, mach)
        stat = run_cmd("{} {} >& funit.log".format(unit_test_tool, args), from_dir=rundir)[0]

        append_testlog(open(os.path.join(rundir, "funit.log"), "r").read())

        expect(stat == 0, "RUN FAIL for FUNIT")

    # Funit is a bit of an oddball test since it's not really running the E3SM model
    # We need to override some methods to make the core infrastructure work.

    def _generate_baseline(self):
        with self._test_status:
            self._test_status.set_status(GENERATE_PHASE, TEST_PASS_STATUS)

    def _compare_baseline(self):
        with self._test_status:
            self._test_status.set_status(BASELINE_PHASE, TEST_PASS_STATUS)
