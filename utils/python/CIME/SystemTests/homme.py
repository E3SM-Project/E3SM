"""
CIME HOMME test. This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.build import post_build

import shutil

logger = logging.getLogger(__name__)

class HOMME(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        if not sharedlib_only:
            # Build HOMME
            srcroot  = self._case.get_value("SRCROOT")
            mach     = self._case.get_value("MACH")
            procs    = self._case.get_value("TOTALPES")
            exeroot  = self._case.get_value("EXEROOT")
            baseline = self._case.get_value("BASELINE_ROOT")
            basegen  = self._case.get_value("BASEGEN_CASE")
            basecmp  = self._case.get_value("BASECMP_CASE")
            generate = self._case.get_value("GENERATE_BASELINE")
            gmake    = self._case.get_value("GMAKE")

            basename = basegen if generate else basecmp
            cmake_cmd = "cmake -C %s/components/homme/cmake/machineFiles/%s.cmake -DUSE_NUM_PROCS=%s %s/components/homme -DHOMME_BASELINE_DIR=%s/%s >& homme.bldlog" % (srcroot, mach, procs, srcroot, baseline, basename)

            run_cmd_no_fail(cmake_cmd, from_dir=exeroot)
            run_cmd_no_fail("%s -j8 >> homme.bldlog 2>&1" % gmake, from_dir=exeroot)

            post_build(self._case, [os.path.join(exeroot, "homme.bldlog")])

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
            run_cmd_no_fail("%s -j 4 baseline >& %s" % (gmake, log), from_dir=exeroot)
            if os.path.isdir(full_baseline_dir):
                shutil.rmtree(full_baseline_dir)
            shutil.copytree(os.path.join(exeroot, "tests", "baseline"), full_baseline_dir)
        elif compare:
            run_cmd_no_fail("%s -j 4 check >& %s" % (gmake, log), from_dir=exeroot)
        else:
            run_cmd_no_fail("%s -j 4 baseline >& %s" % (gmake, log), from_dir=exeroot)
