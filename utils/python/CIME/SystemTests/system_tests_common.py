"""
Base class for CIME system tests
"""
import shutil, glob, gzip, time
from CIME.XML.standard_module_setup import *
from CIME.XML.env_run import EnvRun
from CIME.utils import append_status
from CIME.case_setup import case_setup
from CIME.case_run import case_run
from CIME.case_st_archive import case_st_archive
from CIME.test_status import *

import CIME.build as build

logger = logging.getLogger(__name__)

class SystemTestsCommon(object):

    def __init__(self, case, expected=None):
        """
        initialize a CIME system test object, if the file LockedFiles/env_run.orig.xml
        does not exist copy the current env_run.xml file.  If it does exist restore values
        changed in a previous run of the test.
        """
        self._case = case
        caseroot = case.get_value("CASEROOT")
        self._caseroot = caseroot
        self._orig_caseroot = caseroot
        self._runstatus = None
        self._casebaseid = self._case.get_value("CASEBASEID")
        self._test_status = TestStatus(test_dir=caseroot, test_name=self._casebaseid)

        # Needed for sh scripts
        os.environ["CASEROOT"] = caseroot

        if os.path.isfile(os.path.join(caseroot, "LockedFiles", "env_run.orig.xml")):
            self.compare_env_run(expected=expected)
        elif os.path.isfile(os.path.join(caseroot, "env_run.xml")):
            lockedfiles = os.path.join(caseroot, "LockedFiles")
            try:
                os.stat(lockedfiles)
            except:
                os.mkdir(lockedfiles)
            shutil.copy(os.path.join(caseroot,"env_run.xml"),
                        os.path.join(lockedfiles, "env_run.orig.xml"))

        if self._case.get_value("IS_FIRST_RUN"):
            self._case.set_initial_test_values()

        case_setup(self._case, reset=True, test_mode=True)
        self._case.set_value("TEST",True)
        self._case.flush()

    def build(self, sharedlib_only=False, model_only=False):
        """
        Do NOT override this method, this method is the framework that
        controls the build phase. build_phase is the extension point
        that subclasses should use.
        """
        success = True
        for phase_name, phase_bool in [(SHAREDLIB_BUILD_PHASE, not model_only),
                                       (MODEL_BUILD_PHASE, not sharedlib_only)]:
            if phase_bool:
                with self._test_status:
                    self._test_status.set_status(phase_name, TEST_PENDING_STATUS)

                start_time = time.time()
                try:
                    self.build_phase(sharedlib_only=(phase_name==SHAREDLIB_BUILD_PHASE),
                                     model_only=(phase_name==MODEL_BUILD_PHASE))
                except:
                    success = False
                    logger.warning("Exception during build:\n%s" % (sys.exc_info()[1]))

                time_taken = time.time() - start_time
                with self._test_status:
                    self._test_status.set_status(phase_name, TEST_PASS_STATUS if success else TEST_FAIL_STATUS, comments=("time=%d" % int(time_taken)))

                if not success:
                    break

        return success

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        This is the default build phase implementation, it just does an individual build.
        This is the subclass' extension point if they need to define a custom build
        phase.

        PLEASE THROW EXCEPTION ON FAIL
        """
        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def build_indv(self, sharedlib_only=False, model_only=False):
        """
        Perform an individual build
        """
        build.case_build(self._caseroot, case=self._case,
                         sharedlib_only=sharedlib_only, model_only=model_only)

    def clean_build(self, comps=None):
        build.clean(self._case, cleanlist=comps)

    def run(self):
        """
        Do NOT override this method, this method is the framework that controls
        the run phase. run_phase is the extension point that subclasses should use.
        """
        success = True
        start_time = time.time()
        try:
            expect(self._test_status.get_status(MODEL_BUILD_PHASE) == TEST_PASS_STATUS,
                   "Model was not built!")
            with self._test_status:
                self._test_status.set_status(RUN_PHASE, TEST_PENDING_STATUS)

            self.run_phase()

            if self._case.get_value("GENERATE_BASELINE"):
                self._generate_baseline()

            if self._case.get_value("COMPARE_BASELINE"):
                self._compare_baseline()

            self._check_for_memleak()

        except:
            success = False
            logger.warning("Exception during run:\n%s" % (sys.exc_info()[1]))

        # Always try to report, should NOT throw an exception
        self.report()

        # Writing the run status should be the very last thing due to wait_for_tests
        time_taken = time.time() - start_time
        status = TEST_PASS_STATUS if success else TEST_FAIL_STATUS
        with self._test_status:
            self._test_status.set_status(RUN_PHASE, status, comments=("time=%d" % int(time_taken)))

        # We only return success if every phase, build and later, passed
        return self._test_status.get_overall_test_status(ignore_namelists=True) == TEST_PASS_STATUS

    def run_phase(self):
        """
        This is the default run phase implementation, it just does an individual run.
        This is the subclass' extension point if they need to define a custom run phase.

        PLEASE THROW AN EXCEPTION ON FAIL
        """
        self.run_indv()

    def _set_active_case(self, case):
        """
        Use for tests that have multiple cases
        """
        self._case = case
        self._caseroot = case.get_value("CASEROOT")

    def run_indv(self, suffix="base", coupler_log_path=None, st_archive=False):
        """
        Perform an individual run. Raises an EXCEPTION on fail.
        """
        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        run_type    = self._case.get_value("RUN_TYPE")
        rest_option = self._case.get_value("REST_OPTION")
        rest_n      = self._case.get_value("REST_N")
        infostr     = "doing an %d %s %s test" % (stop_n, stop_option,run_type)

        if rest_option == "none":
            infostr += ", no restarts written"
        else:
            infostr += ", with restarts every %d %s"%(rest_n, rest_option)
        logger.info(infostr)

        case_run(self._case)
        if st_archive:
            case_st_archive(self._case)

        if not self._coupler_log_indicates_run_complete(coupler_log_path):
            expect(False, "Coupler did not indicate run passed")

        if suffix is not None:
            self._component_compare_move(suffix)

    def _coupler_log_indicates_run_complete(self, coupler_log_path):
        newestcpllogfile = self._get_latest_cpl_log(coupler_log_path)
        logger.debug("Latest Coupler log file is %s" % newestcpllogfile)
        # Exception is raised if the file is not compressed
        try:
            if "SUCCESSFUL TERMINATION" in gzip.open(newestcpllogfile, 'rb').read():
                return True
        except:
            logger.info("%s is not compressed, assuming run failed"%newestcpllogfile)
        return False

    def report(self):
        """
        Please explain what kind of things happen in report
        """
        pass

    def _component_compare_move(self, suffix):
        cmd = os.path.join(self._case.get_value("SCRIPTSROOT"), "Tools",
                           "component_compare_move.sh")
        rc, out, err = run_cmd("%s -rundir %s -testcase %s -suffix %s" %
                               (cmd, self._case.get_value('RUNDIR'), self._case.get_value('CASE'), suffix))
        if rc == 0:
            append_status(out, sfile="TestStatus.log")
        else:
            append_status("Component_compare_move.sh failed out: %s\n\nerr: %s\n" % (out, err),
                          sfile="TestStatus.log")

    def _component_compare_test(self, suffix1, suffix2):
        """
        Return value is not generally checked, but is provided in case a custom
        run case needs indirection based on success.
        """
        cmd = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools",
                           "component_compare_test.sh")
        rc, out, err = run_cmd("%s -rundir %s -testcase %s -testcase_base %s -suffix1 %s -suffix2 %s -msg 'Compare %s and %s'"
                               %(cmd, self._case.get_value('RUNDIR'), self._case.get_value('CASE'),
                                 self._case.get_value('CASEBASEID'), suffix1, suffix2, suffix1, suffix2))
        logger.debug("run %s results %d %s %s"%(cmd,rc,out,err))
        status = TEST_PASS_STATUS if rc == 0 else TEST_FAIL_STATUS
        with self._test_status:
            self._test_status.set_status("%s_%s_%s" % (COMPARE_PHASE, suffix1, suffix2), status)

        if rc != 0:
            append_status("Component_compare_test.sh failed out: %s\n\nerr: %s\n"%(out,err),
                          sfile="TestStatus.log")
            return False

        return True

    def _get_mem_usage(self, cpllog):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """
        memlist = []
        meminfo = re.compile(r".*model date =\s+(\w+).*memory =\s+(\d+\.?\d+).*highwater")
        if cpllog is not None and os.path.isfile(cpllog):
            with gzip.open(cpllog, "rb") as f:
                for line in f:
                    m = meminfo.match(line)
                    if m:
                        memlist.append((float(m.group(1)), float(m.group(2))))
        return memlist

    def _get_throughput(self, cpllog):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """
        if cpllog is not None and os.path.isfile(cpllog):
            with gzip.open(cpllog, "rb") as f:
                cpltext = f.read()
                m = re.search(r"# simulated years / cmp-day =\s+(\d+\.\d+)\s",cpltext)
                if m:
                    return float(m.group(1))
        return None

    def _check_for_memleak(self):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """
        cpllog = self._get_latest_cpl_log()

        memlist = self._get_mem_usage(cpllog)

        with self._test_status:
            if len(memlist)<3:
                self._test_status.set_status(MEMLEAK_PHASE, TEST_PASS_STATUS, comments="insuffiencient data for memleak test")
            else:
                finaldate = int(memlist[-1][0])
                originaldate = int(memlist[0][0])
                finalmem = float(memlist[-1][1])
                originalmem = float(memlist[0][1])
                memdiff = -1
                if originalmem > 0:
                    memdiff = (finalmem - originalmem)/originalmem

                if memdiff < 0:
                    self._test_status.set_status(MEMLEAK_PHASE, TEST_PASS_STATUS, comments="insuffiencient data for memleak test")
                elif memdiff < 0.1:
                    self._test_status.set_status(MEMLEAK_PHASE, TEST_PASS_STATUS)
                else:
                    comment = "memleak detected, memory went from %f to %f in %d days" % (originalmem, finalmem, finaldate-originaldate)
                    append_status(comment, sfile="TestStatus.log")
                    self._test_status.set_status(MEMLEAK_PHASE, TEST_FAIL_STATUS, comments=comment)

    def compare_env_run(self, expected=None):
        f1obj = EnvRun(self._caseroot, "env_run.xml")
        f2obj = EnvRun(self._caseroot, os.path.join("LockedFiles", "env_run.orig.xml"))
        diffs = f1obj.compare_xml(f2obj)
        for key in diffs.keys():
            if expected is not None and key in expected:
                logging.warn("  Resetting %s for test"%key)
                f1obj.set_value(key, f2obj.get_value(key, resolved=False))
            else:
                print "Found difference in %s: case: %s original value %s" %\
                    (key, diffs[key][0], diffs[key][1])
                print " Use option --force to run the test with this"\
                    " value or --reset to reset to original"
                return False
        return True

    def _get_latest_cpl_log(self, coupler_log_path=None):
        """
        find and return the latest cpl log file in the run directory
        """
        coupler_log_path = self._case.get_value("RUNDIR") if coupler_log_path is None else coupler_log_path
        cpllog = None
        cpllogs = glob.glob(os.path.join(coupler_log_path, 'cpl.log.*'))
        if cpllogs:
            cpllog = max(cpllogs, key=os.path.getctime)

        return cpllog

    def _compare_baseline(self):
        """
        compare the current test output to a baseline result
        """
        with self._test_status:
            baselineroot = self._case.get_value("BASELINE_ROOT")
            basecmp_dir = os.path.join(baselineroot, self._case.get_value("BASECMP_CASE"))
            for bdir in (baselineroot, basecmp_dir):
                if not os.path.isdir(bdir):
                    comment = "ERROR %s does not exist" % bdir
                    self._test_status.set_status("%s_baseline" % COMPARE_PHASE, TEST_FAIL_STATUS, comments=comment)
                    append_status(comment, sfile="TestStatus.log")
                    return -1

            compgen = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools",
                                   "component_compgen_baseline.sh")
            compgen += " -baseline_dir "+basecmp_dir
            compgen += " -test_dir "+self._case.get_value("RUNDIR")
            compgen += " -compare_tag "+self._case.get_value("BASELINE_NAME_CMP")
            compgen += " -testcase "+self._case.get_value("CASE")
            compgen += " -testcase_base "+self._case.get_value("CASEBASEID")
            rc, out, err = run_cmd(compgen)

            status = TEST_PASS_STATUS if rc == 0 else TEST_FAIL_STATUS
            self._test_status.set_status("%s_baseline" % COMPARE_PHASE, status)
            append_status("Baseline compare results: %s\n%s"%(out,err), sfile="TestStatus.log")

            # compare memory usage to baseline
            newestcpllogfile = self._get_latest_cpl_log()
            memlist = self._get_mem_usage(newestcpllogfile)
            baselog = os.path.join(basecmp_dir, "cpl.log.gz")
            if not os.path.isfile(baselog):
                # for backward compatibility
                baselog = os.path.join(basecmp_dir, "cpl.log")
            if os.path.isfile(baselog) and len(memlist) > 3:
                blmem = self._get_mem_usage(baselog)[-1][1]
                curmem = memlist[-1][1]
                diff = (curmem-blmem)/blmem
                if(diff < 0.1):
                    self._test_status.set_status(MEMCOMP_PHASE, TEST_PASS_STATUS)
                else:
                    comment = "Error: Memory usage increase > 10% from baseline"
                    self._test_status.set_status(MEMCOMP_PHASE, TEST_FAIL_STATUS, comments=comment)
                    append_status(comment, sfile="TestStatus.log")

            # compare throughput to baseline
            current = self._get_throughput(newestcpllogfile)
            baseline = self._get_throughput(baselog)
            #comparing ypd so bigger is better
            if baseline is not None and current is not None:
                diff = (baseline - current)/baseline
                if(diff < 0.25):
                    self._test_status.set_status(THROUGHPUT_PHASE, TEST_PASS_STATUS)
                else:
                    comment = "Error: Computation time increase > 25% from baseline"
                    self._test_status.set_status(THROUGHPUT_PHASE, TEST_FAIL_STATUS, comments=comment)
                    append_status(comment, sfile="TestStatus.log")

    def _generate_baseline(self):
        """
        generate a new baseline case based on the current test
        """
        with self._test_status:
            newestcpllogfile = self._get_latest_cpl_log()
            baselineroot = self._case.get_value("BASELINE_ROOT")
            basegen_dir = os.path.join(baselineroot, self._case.get_value("BASEGEN_CASE"))
            for bdir in (baselineroot, basegen_dir):
                if not os.path.isdir(bdir):
                    comment = "ERROR %s does not exist" % bdir
                    self._test_status.set_status("%s" % GENERATE_PHASE, TEST_FAIL_STATUS, comments=comment)
                    append_status(comment, sfile="TestStatus.log")
                    return -1

            compgen = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools",
                                   "component_compgen_baseline.sh")
            compgen += " -baseline_dir "+basegen_dir
            compgen += " -test_dir "+self._case.get_value("RUNDIR")
            compgen += " -generate_tag "+self._case.get_value("BASELINE_NAME_GEN")
            compgen += " -testcase "+self._case.get_value("CASE")
            compgen += " -testcase_base "+self._case.get_value("CASEBASEID")
            rc, out, err = run_cmd(compgen)

            status = TEST_PASS_STATUS if rc == 0 else TEST_FAIL_STATUS
            self._test_status.set_status("%s" % GENERATE_PHASE, status)
            append_status("Baseline generate results: %s\n%s"%(out,err), sfile="TestStatus.log")

            # copy latest cpl log to baseline
            # drop the date so that the name is generic
            shutil.copyfile(newestcpllogfile,
                            os.path.join(basegen_dir,"cpl.log.gz"))

class FakeTest(SystemTestsCommon):
    """
    Inheriters of the FakeTest Class are intended to test the code.

    All members of the FakeTest Class must
    have names beginnig with "TEST" this is so that the find_system_test
    in utils.py will work with these classes.
    """
    def _set_script(self, script):
        self._script = script # pylint: disable=attribute-defined-outside-init

    def build_phase(self, sharedlib_only=False, model_only=False):
        if (not sharedlib_only):
            exeroot = self._case.get_value("EXEROOT")
            cime_model = self._case.get_value("MODEL")
            modelexe = os.path.join(exeroot, "%s.exe" % cime_model)

            with open(modelexe, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(self._script)

            os.chmod(modelexe, 0755)

            build.post_build(self._case, [])

    def run_phase(self):
        self.run_indv(suffix=None)

class TESTRUNPASS(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        script = \
"""
echo Insta pass
echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID
""" % rundir
        self._set_script(script)
        FakeTest.build_phase(self,
                       sharedlib_only=sharedlib_only, model_only=model_only)

class TESTRUNDIFF(FakeTest):
    """
    You can generate a diff with this test as follows:
    1) Run the test and generate a baseline
    2) set TESTRUNDIFF_ALTERNATE environment variable to TRUE
    3) Re-run the same test from step 1 but do a baseline comparison instead of generation
      3.a) This should give you a DIFF
    """
    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        script = \
"""
echo Insta pass
echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID
if [ -z "$TESTRUNDIFF_ALTERNATE" ]; then
  cp %s/utils/python/tests/cpl.hi1.nc.test %s/%s.cpl.hi.0.nc.base
else
  cp %s/utils/python/tests/cpl.hi2.nc.test %s/%s.cpl.hi.0.nc.base
fi
""" % (rundir, cimeroot, rundir, case, cimeroot, rundir, case)
        self._set_script(script)
        FakeTest.build_phase(self,
                       sharedlib_only=sharedlib_only, model_only=model_only)

class TESTRUNFAIL(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        script = \
"""
echo Insta fail
echo model failed > %s/cpl.log.$LID
exit -1
""" % rundir
        self._set_script(script)
        FakeTest.build_phase(self,
                             sharedlib_only=sharedlib_only, model_only=model_only)

class TESTBUILDFAIL(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        if (not sharedlib_only):
            expect(False, "ERROR: Intentional fail for testing infrastructure")

class TESTRUNSLOWPASS(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        script = \
"""
sleep 300
echo Slow pass
echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID
""" % rundir
        self._set_script(script)
        FakeTest.build_phase(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)

class TESTMEMLEAKFAIL(FakeTest):
    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        testfile = os.path.join(cimeroot,"utils","python","tests","cpl.log.failmemleak.gz")
        script = \
"""
echo Insta pass
gunzip -c %s > %s/cpl.log.$LID
""" % (testfile, rundir)
        self._set_script(script)
        FakeTest.build_phase(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)

class TESTMEMLEAKPASS(FakeTest):
    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        testfile = os.path.join(cimeroot,"utils","python","tests","cpl.log.passmemleak.gz")
        script = \
"""
echo Insta pass
gunzip -c %s > %s/cpl.log.$LID
""" % (testfile, rundir)
        self._set_script(script)
        FakeTest.build_phase(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)
