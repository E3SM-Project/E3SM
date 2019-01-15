"""
Base class for CIME system tests
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.env_run import EnvRun
from CIME.utils import append_testlog, get_model, safe_copy, get_timestamp
from CIME.test_status import *
from CIME.hist_utils import *
from CIME.provenance import save_test_time
from CIME.locked_files import LOCKED_DIR, lock_file, is_locked
import CIME.build as build

import glob, gzip, time, traceback, six

logger = logging.getLogger(__name__)

class SystemTestsCommon(object):

    def __init__(self, case, expected=None):
        """
        initialize a CIME system test object, if the locked env_run.orig.xml
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
        self._init_environment(caseroot)
        self._init_locked_files(caseroot, expected)
        self._skip_pnl = False
        self._cpllog = "med" if self._case.get_value("COMP_INTERFACE")=="nuopc" else "cpl"

    def _init_environment(self, caseroot):
        """
        Do initializations of environment variables that are needed in __init__
        """
        # Needed for sh scripts
        os.environ["CASEROOT"] = caseroot

    def _init_locked_files(self, caseroot, expected):
        """
        If the locked env_run.orig.xml does not exist, copy the current
        env_run.xml file. If it does exist, restore values changed in a previous
        run of the test.
        """
        if is_locked("env_run.orig.xml"):
            self.compare_env_run(expected=expected)
        elif os.path.isfile(os.path.join(caseroot, "env_run.xml")):
            lock_file("env_run.xml", caseroot=caseroot, newname="env_run.orig.xml")

    def _resetup_case(self, phase, reset=False):
        """
        Re-setup this case. This is necessary if user is re-running an already-run
        phase.
        """
        # We never want to re-setup if we're doing the resubmitted run
        phase_status = self._test_status.get_status(phase)
        if reset or (self._case.get_value("IS_FIRST_RUN") and phase_status != TEST_PEND_STATUS):

            logging.warning("Resetting case due to detected re-run of phase {}".format(phase))
            self._case.set_initial_test_values()

            self._case.case_setup(reset=True, test_mode=True)

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
                self._resetup_case(phase_name)
                with self._test_status:
                    self._test_status.set_status(phase_name, TEST_PEND_STATUS)

                start_time = time.time()
                try:
                    self.build_phase(sharedlib_only=(phase_name==SHAREDLIB_BUILD_PHASE),
                                     model_only=(phase_name==MODEL_BUILD_PHASE))
                except BaseException as e:
                    success = False
                    msg = e.__str__()
                    if "FAILED, cat" in msg or "BUILD FAIL" in msg:
                        # Don't want to print stacktrace for a build failure since that
                        # is not a CIME/infrastructure problem.
                        excmsg = msg
                    else:
                        excmsg = "Exception during build:\n{}\n{}".format(msg, traceback.format_exc())

                    logger.warning(excmsg)
                    append_testlog(excmsg)

                time_taken = time.time() - start_time
                with self._test_status:
                    self._test_status.set_status(phase_name, TEST_PASS_STATUS if success else TEST_FAIL_STATUS, comments=("time={:d}".format(int(time_taken))))

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
        model = self._case.get_value('MODEL')
        build.case_build(self._caseroot, case=self._case,
                         sharedlib_only=sharedlib_only, model_only=model_only,
                         save_build_provenance=not model=='cesm')

    def clean_build(self, comps=None):
        if comps is None:
            comps = [x.lower() for x in self._case.get_values("COMP_CLASSES")]
        build.clean(self._case, cleanlist=comps)

    def run(self, skip_pnl=False):
        """
        Do NOT override this method, this method is the framework that controls
        the run phase. run_phase is the extension point that subclasses should use.
        """
        success = True
        start_time = time.time()
        self._skip_pnl = skip_pnl
        try:
            self._resetup_case(RUN_PHASE)
            with self._test_status:
                self._test_status.set_status(RUN_PHASE, TEST_PEND_STATUS)

            self.run_phase()

            if self._case.get_value("GENERATE_BASELINE"):
                self._phase_modifying_call(GENERATE_PHASE, self._generate_baseline)

            if self._case.get_value("COMPARE_BASELINE"):
                self._phase_modifying_call(BASELINE_PHASE,   self._compare_baseline)
                self._phase_modifying_call(MEMCOMP_PHASE,    self._compare_memory)
                self._phase_modifying_call(THROUGHPUT_PHASE, self._compare_throughput)

            self._phase_modifying_call(MEMLEAK_PHASE, self._check_for_memleak)

            self._phase_modifying_call(STARCHIVE_PHASE, self._st_archive_case_test)

        except BaseException as e:
            success = False
            msg = e.__str__()
            if "RUN FAIL" in msg:
                # Don't want to print stacktrace for a model failure since that
                # is not a CIME/infrastructure problem.
                excmsg = msg
            else:
                excmsg = "Exception during run:\n{}\n{}".format(msg, traceback.format_exc())
            logger.warning(excmsg)
            append_testlog(excmsg)

        # Writing the run status should be the very last thing due to wait_for_tests
        time_taken = time.time() - start_time
        status = TEST_PASS_STATUS if success else TEST_FAIL_STATUS
        with self._test_status:
            self._test_status.set_status(RUN_PHASE, status, comments=("time={:d}".format(int(time_taken))))

        if success and get_model() == "e3sm":
            save_test_time(self._case.get_value("BASELINE_ROOT"), self._casebaseid, time_taken)

        if get_model() == "cesm" and self._case.get_value("GENERATE_BASELINE"):
            baseline_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), self._case.get_value("BASEGEN_CASE"))
            generate_teststatus(self._caseroot, baseline_dir)

        # We return success if the run phase worked; memleaks, diffs will not be taken into account
        # with this return value.
        return success

    def run_phase(self):
        """
        This is the default run phase implementation, it just does an individual run.
        This is the subclass' extension point if they need to define a custom run phase.

        PLEASE THROW AN EXCEPTION ON FAIL
        """
        self.run_indv()

    def _get_caseroot(self):
        """
        Returns the current CASEROOT value
        """
        return self._caseroot

    def _set_active_case(self, case):
        """
        Use for tests that have multiple cases
        """
        self._case = case
        self._case.load_env(reset=True)
        self._caseroot = case.get_value("CASEROOT")

    def run_indv(self, suffix="base", st_archive=False):
        """
        Perform an individual run. Raises an EXCEPTION on fail.
        """
        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        run_type    = self._case.get_value("RUN_TYPE")
        rundir      = self._case.get_value("RUNDIR")
        is_batch    = self._case.get_value("BATCH_SYSTEM") != "none"

        # remove any cprnc output leftover from previous runs
        for compout in glob.iglob(os.path.join(rundir,"*.cprnc.out")):
            os.remove(compout)

        infostr     = "doing an {:d} {} {} test".format(stop_n, stop_option, run_type)

        rest_option = self._case.get_value("REST_OPTION")
        if rest_option == "none" or rest_option == "never":
            infostr += ", no restarts written"
        else:
            rest_n   = self._case.get_value("REST_N")
            infostr += ", with restarts every {:d} {}".format(rest_n, rest_option)

        logger.info(infostr)

        self._case.case_run(skip_pnl=self._skip_pnl, submit_resubmits=is_batch)

        if not self._coupler_log_indicates_run_complete():
            expect(False, "Coupler did not indicate run passed")

        if suffix is not None:
            self._component_compare_copy(suffix)

        if st_archive:
            self._case.case_st_archive(resubmit=True)

    def _coupler_log_indicates_run_complete(self):
        newestcpllogfiles = self._get_latest_cpl_logs()
        logger.debug("Latest Coupler log file(s) {}" .format(newestcpllogfiles))
        # Exception is raised if the file is not compressed
        allgood = len(newestcpllogfiles)
        for cpllog in newestcpllogfiles:
            try:
                if six.b("SUCCESSFUL TERMINATION") in gzip.open(cpllog, 'rb').read():
                    allgood = allgood - 1
            except BaseException as e:
                msg = e.__str__()

                logger.info("{} is not compressed, assuming run failed {}".format(cpllog, msg))

        return allgood==0

    def _component_compare_copy(self, suffix):
        comments = copy(self._case, suffix)
        append_testlog(comments)

    def _component_compare_test(self, suffix1, suffix2, success_change=False):
        """
        Return value is not generally checked, but is provided in case a custom
        run case needs indirection based on success.
        If success_change is True, success requires some files to be different
        """
        success, comments = self._do_compare_test(suffix1, suffix2)
        if success_change:
            success = not success

        append_testlog(comments)
        status = TEST_PASS_STATUS if success else TEST_FAIL_STATUS
        with self._test_status:
            self._test_status.set_status("{}_{}_{}".format(COMPARE_PHASE, suffix1, suffix2), status)
        return success

    def _do_compare_test(self, suffix1, suffix2):
        """
        Wraps the call to compare_test to facilitate replacement in unit
        tests
        """
        return compare_test(self._case, suffix1, suffix2)

    def _st_archive_case_test(self):
        result = self._case.test_env_archive()
        with self._test_status:
            if result:
                self._test_status.set_status(STARCHIVE_PHASE, TEST_PASS_STATUS)
            else:
                self._test_status.set_status(STARCHIVE_PHASE, TEST_FAIL_STATUS)

    def _get_mem_usage(self, cpllog):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """
        memlist = []
        meminfo = re.compile(r".*model date =\s+(\w+).*memory =\s+(\d+\.?\d+).*highwater")
        if cpllog is not None and os.path.isfile(cpllog):
            if '.gz' == cpllog[-3:]:
                fopen = gzip.open
            else:
                fopen = open
            with fopen(cpllog, "rb") as f:
                for line in f:
                    m = meminfo.match(line.decode('utf-8'))
                    if m:
                        memlist.append((float(m.group(1)), float(m.group(2))))
        # Remove the last mem record, it's sometimes artificially high
        if len(memlist) > 0:
            memlist.pop()
        return memlist

    def _get_throughput(self, cpllog):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """
        if cpllog is not None and os.path.isfile(cpllog):
            with gzip.open(cpllog, "rb") as f:
                cpltext = f.read().decode('utf-8')
                m = re.search(r"# simulated years / cmp-day =\s+(\d+\.\d+)\s",cpltext)
                if m:
                    return float(m.group(1))
        return None

    def _phase_modifying_call(self, phase, function):
        """
        Ensures that unexpected exceptions from phases will result in a FAIL result
        in the TestStatus file for that phase.
        """
        try:
            function()
        except BaseException as e:
            msg = e.__str__()
            excmsg = "Exception during {}:\n{}\n{}".format(phase, msg, traceback.format_exc())

            logger.warning(excmsg)
            append_testlog(excmsg)

            with self._test_status:
                self._test_status.set_status(phase, TEST_FAIL_STATUS, comments="exception")

    def _check_for_memleak(self):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """
        with self._test_status:
            latestcpllogs = self._get_latest_cpl_logs()
            for cpllog in latestcpllogs:
                memlist = self._get_mem_usage(cpllog)

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
                    tolerance = self._case.get_value("TEST_MEMLEAK_TOLERANCE")
                    if tolerance is None:
                        tolerance = 0.1
                    expect(tolerance > 0.0, "Bad value for memleak tolerance in test")
                    if memdiff < 0:
                        self._test_status.set_status(MEMLEAK_PHASE, TEST_PASS_STATUS, comments="insuffiencient data for memleak test")
                    elif memdiff < tolerance:
                        self._test_status.set_status(MEMLEAK_PHASE, TEST_PASS_STATUS)
                    else:
                        comment = "memleak detected, memory went from {:f} to {:f} in {:d} days".format(originalmem, finalmem, finaldate-originaldate)
                        append_testlog(comment)
                        self._test_status.set_status(MEMLEAK_PHASE, TEST_FAIL_STATUS, comments=comment)

    def compare_env_run(self, expected=None):
        """
        Compare env_run file to original and warn about differences
        """
        components = self._case.get_values("COMP_CLASSES")
        f1obj = self._case.get_env("run")
        f2obj = EnvRun(self._caseroot, os.path.join(LOCKED_DIR, "env_run.orig.xml"), components=components)
        diffs = f1obj.compare_xml(f2obj)
        for key in diffs.keys():
            if expected is not None and key in expected:
                logging.warning("  Resetting {} for test".format(key))
                f1obj.set_value(key, f2obj.get_value(key, resolved=False))
            else:
                print("WARNING: Found difference in test {}: case: {} original value {}".format(key, diffs[key][0], diffs[key][1]))
                return False
        return True

    def _get_latest_cpl_logs(self):
        """
        find and return the latest cpl log file in the run directory
        """
        coupler_log_path = self._case.get_value("RUNDIR")
        cpllogs = glob.glob(os.path.join(coupler_log_path, '{}*.log.*'.format(self._cpllog)))
        lastcpllogs = []
        if cpllogs:
            lastcpllogs.append(max(cpllogs, key=os.path.getctime))
            basename = os.path.basename(lastcpllogs[0])
            suffix = basename.split('.',1)[1]
            for log in cpllogs:
                if log in lastcpllogs:
                    continue

                if log.endswith(suffix):
                    lastcpllogs.append(log)

        return lastcpllogs

    def _compare_memory(self):
        with self._test_status:
            # compare memory usage to baseline
            baseline_name = self._case.get_value("BASECMP_CASE")
            basecmp_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), baseline_name)
            newestcpllogfiles = self._get_latest_cpl_logs()
            if len(newestcpllogfiles) > 0:
                memlist = self._get_mem_usage(newestcpllogfiles[0])
            for cpllog in newestcpllogfiles:
                m = re.search(r"/({}.*.log).*.gz".format(self._cpllog),cpllog)
                if m is not None:
                    baselog = os.path.join(basecmp_dir, m.group(1))+".gz"
                if baselog is None or not os.path.isfile(baselog):
                    # for backward compatibility
                    baselog = os.path.join(basecmp_dir, self._cpllog+".log")
                if os.path.isfile(baselog) and len(memlist) > 3:
                    blmem = self._get_mem_usage(baselog)
                    blmem = 0 if blmem == [] else blmem[-1][1]
                    curmem = memlist[-1][1]
                    diff = (curmem-blmem)/blmem
                    if diff < 0.1 and self._test_status.get_status(MEMCOMP_PHASE) is None:
                        self._test_status.set_status(MEMCOMP_PHASE, TEST_PASS_STATUS)
                    elif self._test_status.get_status(MEMCOMP_PHASE) != TEST_FAIL_STATUS:
                        comment = "Error: Memory usage increase > 10% from baseline"
                        self._test_status.set_status(MEMCOMP_PHASE, TEST_FAIL_STATUS, comments=comment)
                        append_testlog(comment)

    def _compare_throughput(self):
        with self._test_status:
            # compare memory usage to baseline
            baseline_name = self._case.get_value("BASECMP_CASE")
            basecmp_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), baseline_name)
            newestcpllogfiles = self._get_latest_cpl_logs()
            for cpllog in newestcpllogfiles:
                m = re.search(r"/({}.*.log).*.gz".format(self._cpllog), cpllog)
                if m is not None:
                    baselog = os.path.join(basecmp_dir, m.group(1))+".gz"
                if baselog is None or not os.path.isfile(baselog):
                    # for backward compatibility
                    baselog = os.path.join(basecmp_dir, self._cpllog)

                if os.path.isfile(baselog):
                    # compare throughput to baseline
                    current = self._get_throughput(cpllog)
                    baseline = self._get_throughput(baselog)
                    #comparing ypd so bigger is better
                    if baseline is not None and current is not None:
                        diff = (baseline - current)/baseline
                        tolerance = self._case.get_value("TEST_TPUT_TOLERANCE")
                        if tolerance is None:
                            tolerance = 0.25
                        expect(tolerance > 0.0, "Bad value for throughput tolerance in test")
                        if diff < tolerance and self._test_status.get_status(THROUGHPUT_PHASE) is None:
                            self._test_status.set_status(THROUGHPUT_PHASE, TEST_PASS_STATUS)
                        elif self._test_status.get_status(THROUGHPUT_PHASE) != TEST_FAIL_STATUS:
                            comment = "Error: Computation time increase > {:d} pct from baseline".format(int(tolerance*100))
                            self._test_status.set_status(THROUGHPUT_PHASE, TEST_FAIL_STATUS, comments=comment)
                            append_testlog(comment)

    def _compare_baseline(self):
        """
        compare the current test output to a baseline result
        """
        with self._test_status:
            # compare baseline
            success, comments = compare_baseline(self._case)
            append_testlog(comments)
            status = TEST_PASS_STATUS if success else TEST_FAIL_STATUS
            baseline_name = self._case.get_value("BASECMP_CASE")
            ts_comments = os.path.dirname(baseline_name) + ": " + get_ts_synopsis(comments)
            self._test_status.set_status(BASELINE_PHASE, status, comments=ts_comments)

    def _generate_baseline(self):
        """
        generate a new baseline case based on the current test
        """
        with self._test_status:
            # generate baseline
            success, comments = generate_baseline(self._case)
            append_testlog(comments)
            status = TEST_PASS_STATUS if success else TEST_FAIL_STATUS
            baseline_name = self._case.get_value("BASEGEN_CASE")
            self._test_status.set_status(GENERATE_PHASE, status, comments=os.path.dirname(baseline_name))
            basegen_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), self._case.get_value("BASEGEN_CASE"))
            # copy latest cpl log to baseline
            # drop the date so that the name is generic
            newestcpllogfiles = self._get_latest_cpl_logs()
            for cpllog in newestcpllogfiles:
                m = re.search(r"/({}.*.log).*.gz".format(self._cpllog),cpllog)
                if m is not None:
                    baselog = os.path.join(basegen_dir, m.group(1))+".gz"
                    safe_copy(cpllog,
                              os.path.join(basegen_dir,baselog))

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
            modelexe = os.path.join(exeroot, "{}.exe".format(cime_model))

            with open(modelexe, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(self._script)

            os.chmod(modelexe, 0o755)

            build.post_build(self._case, [], build_complete=True)

    def run_indv(self, suffix="base", st_archive=False):
        mpilib = self._case.get_value("MPILIB")
        # This flag is needed by mpt to run a script under mpiexec
        if mpilib == "mpt":
            os.environ["MPI_SHEPHERD"] = "true"
        super(FakeTest, self).run_indv(suffix, st_archive)

class TESTRUNPASS(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        script = \
"""
echo Insta pass
echo SUCCESSFUL TERMINATION > {rundir}/{log}.log.$LID
cp {root}/scripts/tests/cpl.hi1.nc.test {rundir}/{case}.cpl.hi.0.nc
""".format(rundir=rundir, log=self._cpllog, root=cimeroot, case=case)
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
echo SUCCESSFUL TERMINATION > {rundir}/{log}.log.$LID
if [ -z "$TESTRUNDIFF_ALTERNATE" ]; then
  cp {root}/scripts/tests/cpl.hi1.nc.test {rundir}/{case}.cpl.hi.0.nc
else
  cp {root}/scripts/tests/cpl.hi2.nc.test {rundir}/{case}.cpl.hi.0.nc
fi
""".format(rundir=rundir, log=self._cpllog, root=cimeroot, case=case)
        self._set_script(script)
        FakeTest.build_phase(self,
                       sharedlib_only=sharedlib_only, model_only=model_only)

class TESTTESTDIFF(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        script = \
"""
echo Insta pass
echo SUCCESSFUL TERMINATION > {rundir}/{log}.log.$LID
cp {root}/scripts/tests/cpl.hi1.nc.test {rundir}/{case}.cpl.hi.0.nc
cp {root}/scripts/tests/cpl.hi2.nc.test {rundir}/{case}.cpl.hi.0.nc.rest
""".format(rundir=rundir, log=self._cpllog, root=cimeroot, case=case)
        self._set_script(script)
        super(TESTTESTDIFF, self).build_phase(sharedlib_only=sharedlib_only,
                                              model_only=model_only)

    def run_phase(self):
        super(TESTTESTDIFF, self).run_phase()
        self._component_compare_test("base", "rest")

class TESTRUNFAIL(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        script = \
"""
if [ -z "$TESTRUNFAIL_PASS" ]; then
  echo Insta fail
  echo model failed > {rundir}/{log}.log.$LID
  exit -1
else
  echo Insta pass
  echo SUCCESSFUL TERMINATION > {rundir}/{log}.log.$LID
  cp {root}/scripts/tests/cpl.hi1.nc.test {rundir}/{case}.cpl.hi.0.nc
fi
""".format(rundir=rundir, log=self._cpllog, root=cimeroot, case=case)
        self._set_script(script)
        FakeTest.build_phase(self,
                             sharedlib_only=sharedlib_only, model_only=model_only)

class TESTRUNFAILEXC(TESTRUNPASS):

    def run_phase(self):
        raise RuntimeError("Exception from run_phase")

class TESTRUNSTARCFAIL(TESTRUNPASS):

    def _st_archive_case_test(self):
        raise RuntimeError("Exception from st archive")

class TESTBUILDFAIL(TESTRUNPASS):

    def build_phase(self, sharedlib_only=False, model_only=False):
        if "TESTBUILDFAIL_PASS" in os.environ:
            TESTRUNPASS.build_phase(self, sharedlib_only, model_only)
        else:
            if (not sharedlib_only):
                blddir = self._case.get_value("EXEROOT")
                bldlog = os.path.join(blddir, "{}.bldlog.{}".format(get_model(), get_timestamp("%y%m%d-%H%M%S")))
                with open(bldlog, "w") as fd:
                    fd.write("BUILD FAIL: Intentional fail for testing infrastructure")

                expect(False, "BUILD FAIL: Intentional fail for testing infrastructure")

class TESTBUILDFAILEXC(FakeTest):

    def __init__(self, case):
        FakeTest.__init__(self, case)
        raise RuntimeError("Exception from init")

class TESTRUNSLOWPASS(FakeTest):

    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        script = \
"""
sleep 300
echo Slow pass
echo SUCCESSFUL TERMINATION > {rundir}/{log}.log.$LID
cp {root}/scripts/tests/cpl.hi1.nc.test {rundir}/{case}.cpl.hi.0.nc
""".format(rundir=rundir, log=self._cpllog, root=cimeroot, case=case)
        self._set_script(script)
        FakeTest.build_phase(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)

class TESTMEMLEAKFAIL(FakeTest):
    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        testfile = os.path.join(cimeroot,"scripts","tests","cpl.log.failmemleak.gz")
        script = \
"""
echo Insta pass
gunzip -c {testfile} > {rundir}/{log}.log.$LID
cp {root}/scripts/tests/cpl.hi1.nc.test {rundir}/{case}.cpl.hi.0.nc
""".format(testfile=testfile, rundir=rundir, log=self._cpllog, root=cimeroot, case=case)
        self._set_script(script)
        FakeTest.build_phase(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)

class TESTMEMLEAKPASS(FakeTest):
    def build_phase(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        testfile = os.path.join(cimeroot,"scripts","tests","cpl.log.passmemleak.gz")
        script = \
"""
echo Insta pass
gunzip -c {testfile} > {rundir}/{log}.log.$LID
cp {root}/scripts/tests/cpl.hi1.nc.test {rundir}/{case}.cpl.hi.0.nc
""".format(testfile=testfile, rundir=rundir, log=self._cpllog, root=cimeroot, case=case)
        self._set_script(script)
        FakeTest.build_phase(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)
