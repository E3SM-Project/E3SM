"""
CIME HOMME test. This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.build import post_build
from CIME.utils import append_testlog, SharedArea
from CIME.test_status import *

import shutil
from distutils import dir_util

logger = logging.getLogger(__name__)

class HommeBase(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)
        case.load_env()
        self.csnd = 'not defined'
        self.cmakesuffix = self.csnd

    def build_phase(self, sharedlib_only=False, model_only=False):
        if not sharedlib_only:
            # Build HOMME
            srcroot  = self._case.get_value("SRCROOT")
            mach     = self._case.get_value("MACH")
            procs    = self._case.get_value("TOTALPES")
            exeroot  = self._case.get_value("EXEROOT")
            baseline = self._case.get_value("BASELINE_ROOT")
            basecmp  = self._case.get_value("BASECMP_CASE")
            compare  = self._case.get_value("COMPARE_BASELINE")
            gmake    = self._case.get_value("GMAKE")
            gmake_j  = self._case.get_value("GMAKE_J")
            cprnc    = self._case.get_value("CCSM_CPRNC")

            if compare:
                basename = basecmp
                baselinedir = baseline
            else:
                basename = ""
                baselinedir = exeroot

            expect(self.cmakesuffix != self.csnd,
                "ERROR in hommebaseclass: Must have cmakesuffix set up")

            cmake_cmd = "cmake -C {0}/components/homme/cmake/machineFiles/{1}{6}.cmake -DUSE_NUM_PROCS={2} {0}/components/homme -DHOMME_BASELINE_DIR={3}/{4} -DCPRNC_DIR={5}/..".format(srcroot, mach, procs, baselinedir, basename, cprnc, self.cmakesuffix)

            run_cmd_no_fail(cmake_cmd, arg_stdout="homme.bldlog", combine_output=True, from_dir=exeroot)
            run_cmd_no_fail("{} -j{} VERBOSE=1 test-execs".format(gmake, gmake_j), arg_stdout="homme.bldlog", combine_output=True, from_dir=exeroot)

            post_build(self._case, [os.path.join(exeroot, "homme.bldlog")], build_complete=True)

    def run_phase(self):

        rundir   = self._case.get_value("RUNDIR")
        exeroot  = self._case.get_value("EXEROOT")
        baseline = self._case.get_value("BASELINE_ROOT")
        compare  = self._case.get_value("COMPARE_BASELINE")
        generate = self._case.get_value("GENERATE_BASELINE")
        basegen  = self._case.get_value("BASEGEN_CASE")
        gmake    = self._case.get_value("GMAKE")

        log = os.path.join(rundir, "homme.log")
        if os.path.exists(log):
            os.remove(log)

        if generate:
            full_baseline_dir = os.path.join(baseline, basegen, "tests", "baseline")
            stat = run_cmd("{} -j 4 baseline".format(gmake), arg_stdout=log, combine_output=True, from_dir=exeroot)[0]
            if stat == 0:
                if os.path.isdir(full_baseline_dir):
                    shutil.rmtree(full_baseline_dir)

                with SharedArea():
                    dir_util.copy_tree(os.path.join(exeroot, "tests", "baseline"), full_baseline_dir, preserve_mode=False)

        elif compare:
            stat = run_cmd("{} -j 4 check".format(gmake), arg_stdout=log, combine_output=True, from_dir=exeroot)[0]

        else:
            stat = run_cmd("{} -j 4 baseline".format(gmake), arg_stdout=log, combine_output=True, from_dir=exeroot)[0]
            if stat == 0:
                stat = run_cmd("{} -j 4 check".format(gmake), arg_stdout=log, combine_output=True, from_dir=exeroot)[0]

        # Add homme.log output to TestStatus.log so that it can
        # appear on the dashboard. Otherwise, the TestStatus.log
        # is pretty useless for this test.
        append_testlog(open(log, "r").read())

        expect(stat == 0, "RUN FAIL for HOMME")

    # Homme is a bit of an oddball test since it's not really running the E3SM model
    # We need to override some methods to make the core infrastructure work.

    def _generate_baseline(self):
        with self._test_status:
            self._test_status.set_status(GENERATE_PHASE, TEST_PASS_STATUS)

    def _compare_baseline(self):
        with self._test_status:
            self._test_status.set_status(BASELINE_PHASE, TEST_PASS_STATUS)
