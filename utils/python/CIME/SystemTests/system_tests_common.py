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
        self._start_time = time.time()
        caseroot = case.get_value("CASEROOT")
        self._caseroot = caseroot
        self._orig_caseroot = caseroot
        self._runstatus = None
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

    def fail_test(self):
        self._runstatus = "FAIL"

    def has_failed(self):
        return self._runstatus == "FAIL"

    def pass_test(self):
        if not self.has_failed():
            self._runstatus = "PASS"

    def has_passed(self):
        return self._runstatus == "PASS"

    def build(self, sharedlib_only=False, model_only=False):
        build.case_build(self._caseroot, case=self._case,
                         sharedlib_only=sharedlib_only, model_only=model_only)

    def clean_build(self, comps=None):
        build.clean(self._case, cleanlist=comps)

    def run(self):
        return self._run()

    def _set_active_case(self, case):
        """
        Use for tests that have multiple cases
        """
        self._case = case
        self._caseroot = case.get_value("CASEROOT")

    def _run(self, suffix="base", coupler_log_path=None, st_archive=False):
        try:
            success = case_run(self._case)
            if success and st_archive:
                success = case_st_archive(self._case)

            if success and self._coupler_log_indicates_run_complete(coupler_log_path):
                self.pass_test()
            else:
                success = False
                self.fail_test()

            if success and suffix is not None:
                self._component_compare_move(suffix)
        except:
            # An exception must not prevent the TestStatus file from
            # being marked FAIL
            success = False
            self.fail_test()
            logger.warning("Exception during run: %s" % (sys.exc_info()[1]))

        return success

    def __del__(self):
        if self._runstatus is not None:
            test_status = os.path.join(self._orig_caseroot, "TestStatus")
            with open(test_status, 'r') as f:
                teststatusfile = f.read()
            li = teststatusfile.rsplit('PEND', 1)
            teststatusfile = self._runstatus.join(li)
            total_time = int(time.time() - self._start_time)
            with open(test_status, 'w') as f:
                f.write(teststatusfile)
                f.write("COMMENT TIME %d\n" % total_time)

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
        newestcpllogfile = self._get_latest_cpl_log()
        self._check_for_memleak(newestcpllogfile)

    def _component_compare_move(self, suffix):
        cmd = os.path.join(self._case.get_value("SCRIPTSROOT"), "Tools",
                           "component_compare_move.sh")
        rc, out, err = run_cmd("%s -rundir %s -testcase %s -suffix %s" %
                               (cmd, self._case.get_value('RUNDIR'), self._case.get_value('CASE'), suffix))
        if rc == 0:
            append_status(out, sfile="TestStatus.log")
        else:
            append_status("Component_compare_test.sh failed out: %s\n\nerr: %s\n"%(out,err)
                          ,sfile="TestStatus.log")

    def _component_compare_test(self, suffix1, suffix2):
        cmd = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools",
                           "component_compare_test.sh")
        rc, out, err = run_cmd("%s -rundir %s -testcase %s -testcase_base %s -suffix1 %s -suffix2 %s -msg 'Compare %s and %s'"
                               %(cmd, self._case.get_value('RUNDIR'), self._case.get_value('CASE'),
                                 self._case.get_value('CASEBASEID'), suffix1, suffix2, suffix1, suffix2))
        logger.debug("run %s results %d %s %s"%(cmd,rc,out,err))
        if rc == 0:
            append_status(out.replace("compare","compare functionality", 1) + "\n",
                          sfile="TestStatus")
        else:
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

    def _check_for_memleak(self, cpllog):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """
        if not self.has_passed():
            append_status("Cannot check memory, test did not pass.\n", sfile="TestStatus.log")
            return

        memlist = self._get_mem_usage(cpllog)

        if len(memlist)<3:
            append_status("COMMENT: insuffiencient data for memleak test",sfile="TestStatus")
        else:
            finaldate = int(memlist[-1][0])
            originaldate = int(memlist[0][0])
            finalmem = float(memlist[-1][1])
            originalmem = float(memlist[0][1])
            memdiff = -1
            if originalmem > 0:
                memdiff = (finalmem - originalmem)/originalmem
            if memdiff < 0:
                append_status("COMMENT: insuffiencient data for memleak test",sfile="TestStatus")
            elif memdiff < 0.1:
                append_status("PASS %s memleak"%(self._case.get_value("CASEBASEID")),
                             sfile="TestStatus")
            else:
                append_status("memleak detected, memory went from %f to %f in %d days"
                             %(originalmem, finalmem, finaldate-originaldate),sfile="TestStatus.log")
                append_status("FAIL %s memleak"%(self._case.get_value("CASEBASEID")),
                             sfile="TestStatus")

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

    def compare_baseline(self):
        """
        compare the current test output to a baseline result
        """
        if not self.has_passed():
            append_status("Cannot compare baselines, test did not pass.\n", sfile="TestStatus.log")
            return

        baselineroot = self._case.get_value("BASELINE_ROOT")
        basecmp_dir = os.path.join(baselineroot, self._case.get_value("BASECMP_CASE"))
        for bdir in (baselineroot, basecmp_dir):
            if not os.path.isdir(bdir):
                append_status("FAIL %s compare\n"%self._case.get_value("CASEBASEID"),
                              sfile="TestStatus")
                append_status("ERROR %s does not exist"%bdir, sfile="TestStatus.log")
                return -1
        compgen = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools",
                               "component_compgen_baseline.sh")
        compgen += " -baseline_dir "+basecmp_dir
        compgen += " -test_dir "+self._case.get_value("RUNDIR")
        compgen += " -compare_tag "+self._case.get_value("BASELINE_NAME_CMP")
        compgen += " -testcase "+self._case.get_value("CASE")
        compgen += " -testcase_base "+self._case.get_value("CASEBASEID")
        rc, out, err = run_cmd(compgen)

        append_status(out.replace("compare","compare baseline", 1),sfile="TestStatus")
        if rc != 0:
            append_status("Error in Baseline compare: %s\n%s"%(out,err), sfile="TestStatus.log")

        # compare memory usage to baseline
        newestcpllogfile = self._get_latest_cpl_log()
        memlist = self._get_mem_usage(newestcpllogfile)
        baselog = os.path.join(basecmp_dir, "cpl.log.gz")
        if not os.path.isfile(baselog):
            # for backward compatibility
            baselog = os.path.join(basecmp_dir, "cpl.log")
        if len(memlist) > 3:
            blmem = self._get_mem_usage(baselog)[-1][1]
            curmem = memlist[-1][1]
            diff = (curmem-blmem)/blmem
            if(diff < 0.1):
                append_status("PASS %s memcomp\n"%self._case.get_value("CASEBASEID"),
                              sfile="TestStatus")
            else:
                append_status("FAIL %s memcomp\n"%self._case.get_value("CASEBASEID"),
                              sfile="TestStatus")
                append_status("Error in memory compare: Memory usage increase > 10% from baseline",
                              sfile="TestStatus.log")
        # compare throughput to baseline
        current = self._get_throughput(newestcpllogfile)
        baseline = self._get_throughput(baselog)
        #comparing ypd so bigger is better
        if baseline is not None and current is not None:
            diff = (baseline - current)/baseline
            if(diff < 0.25):
                append_status("PASS %s tputcomp\n"%self._case.get_value("CASEBASEID"),
                              sfile="TestStatus")
            else:
                append_status("FAIL %s tputcomp\n"%self._case.get_value("CASEBASEID"),
                              sfile="TestStatus")
                append_status("Error in throughput compare: Computation time increase > 25% from baseline",
                              sfile="TestStatus.log")

    def generate_baseline(self):
        """
        generate a new baseline case based on the current test
        """
        if not self.has_passed():
            append_status("Cannot generate baselines, test did not pass.\n", sfile="TestStatus.log")
            return

        newestcpllogfile = self._get_latest_cpl_log()
        baselineroot = self._case.get_value("BASELINE_ROOT")
        basegen_dir = os.path.join(baselineroot, self._case.get_value("BASEGEN_CASE"))
        for bdir in (baselineroot, basegen_dir):
            if not os.path.isdir(bdir):
                append_status("FAIL %s generate\n" % self._case.get_value("CASEBASEID"),
                             sfile="TestStatus")
                append_status("ERROR %s does not exist" % bdir, sfile="TestStatus.log")
                return -1
        compgen = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools",
                               "component_compgen_baseline.sh")
        compgen += " -baseline_dir "+basegen_dir
        compgen += " -test_dir "+self._case.get_value("RUNDIR")
        compgen += " -generate_tag "+self._case.get_value("BASELINE_NAME_GEN")
        compgen += " -testcase "+self._case.get_value("CASE")
        compgen += " -testcase_base "+self._case.get_value("CASEBASEID")
        rc, out, err = run_cmd(compgen)
        # copy latest cpl log to baseline
        # drop the date so that the name is generic
        shutil.copyfile(newestcpllogfile,
                        os.path.join(basegen_dir,"cpl.log.gz"))
        append_status(out,sfile="TestStatus")
        if rc != 0:
            append_status("Error in Baseline Generate: %s"%err,sfile="TestStatus.log")

class FakeTest(SystemTestsCommon):

    def _set_script(self, script):
        self._script = script # pylint: disable=attribute-defined-outside-init

    def build(self, sharedlib_only=False, model_only=False):
        if (not sharedlib_only):
            exeroot = self._case.get_value("EXEROOT")
            cime_model = self._case.get_value("MODEL")
            modelexe = os.path.join(exeroot, "%s.exe" % cime_model)

            with open(modelexe, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(self._script)

            os.chmod(modelexe, 0755)
            self._case.set_value("BUILD_COMPLETE", True)
            self._case.flush()

    def run(self):
        return SystemTestsCommon._run(self, suffix=None)

class TESTRUNPASS(FakeTest):

    def build(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        script = \
"""
echo Insta pass
echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID
""" % rundir
        self._set_script(script)
        FakeTest.build(self,
                       sharedlib_only=sharedlib_only, model_only=model_only)

class TESTRUNDIFF(FakeTest):

    def build(self, sharedlib_only=False, model_only=False):
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
        FakeTest.build(self,
                       sharedlib_only=sharedlib_only, model_only=model_only)

class TESTRUNFAIL(FakeTest):

    def build(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        script = \
"""
echo Insta fail
echo model failed > %s/cpl.log.$LID
exit -1
""" % rundir
        self._set_script(script)
        FakeTest.build(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)

class TESTBUILDFAIL(FakeTest):

    def build(self, sharedlib_only=False, model_only=False):
        if (not sharedlib_only):
            expect(False, "ERROR: Intentional fail for testing infrastructure")

class TESTRUNSLOWPASS(FakeTest):

    def build(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        script = \
"""
sleep 300
echo Slow pass
echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID
""" % rundir
        self._set_script(script)
        FakeTest.build(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)

class TESTMEMLEAKFAIL(FakeTest):
    def build(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        testfile = os.path.join(cimeroot,"utils","python","tests","cpl.log.failmemleak.gz")
        script = \
"""
echo Insta pass
gunzip -c %s > %s/cpl.log.$LID
""" % (testfile, rundir)
        self._set_script(script)
        FakeTest.build(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)

class TESTMEMLEAKPASS(FakeTest):
    def build(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        testfile = os.path.join(cimeroot,"utils","python","tests","cpl.log.passmemleak.gz")
        script = \
"""
echo Insta pass
gunzip -c %s > %s/cpl.log.$LID
""" % (testfile, rundir)
        self._set_script(script)
        FakeTest.build(self,
                        sharedlib_only=sharedlib_only, model_only=model_only)
