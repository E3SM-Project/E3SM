"""
Base class for CIME system tests
"""
import shutil, glob
from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.XML.env_run import EnvRun
from CIME.utils import run_cmd
from CIME.case_setup import case_setup
import CIME.build as build

logger = logging.getLogger(__name__)

class SystemTestsCommon(object):
    def __init__(self, caseroot=os.getcwd(), case=None):
        """
        initialize a CIME system test object, if the file LockedFiles/env_run.orig.xml
        does not exist copy the current env_run.xml file.  If it does exist restore values
        changed in a previous run of the test.
        """
        print caseroot
        self._caseroot = caseroot
        # Needed for sh scripts
        os.environ["CASEROOT"] = caseroot
        if case is None:
            self._case = Case(caseroot)
        else:
            self._case = case

        if os.path.isfile(os.path.join(caseroot, "LockedFiles", "env_run.orig.xml")):
            self.compare_env_run()
        elif os.path.isfile(os.path.join(caseroot, "env_run.xml")):
            lockedfiles = os.path.join(caseroot, "Lockedfiles")
            try:
                os.stat(lockedfiles)
            except:
                os.mkdir(lockedfiles)
            shutil.copy("env_run.xml",
                        os.path.join(lockedfiles, "env_run.orig.xml"))

        self._case.set_initial_test_values()
        case_setup(self._caseroot, reset=True, test_mode=True)


    def build(self, sharedlib_only=False, model_only=False):
        build.case_build(self._caseroot, case=self._case,
                         sharedlib_only=sharedlib_only, model_only=model_only)


    def run(self):
        with open("TestStatus", 'a') as f:
            f.write("PEND %s RUN\n"%self._case.get_value("CASEBASEID"))

        with open("TestStatus", 'r') as f:
            teststatusfile = f.read()

        rc, out, err = run_cmd("./case.run", ok_to_fail=True)
        if rc == 0:
            result = "PASS"
        else:
            result = "FAIL"

        with open("TestStatus.log", 'a') as f:
            f.write("case.run output is %s"%out)
            f.write("case.run error is %s"%err)

        li = teststatusfile.rsplit('PEND', 1)
        teststatusfile = result.join(li)
        with open("TestStatus", 'w') as f:
            f.write(teststatusfile)


        return

    def report(self):
        self._checkformemleak()
        return

    def _checkformemleak(self):
        """
        Examine memory usage as recorded in the cpl log file and look for unexpected
        increases.
        """

        newestcpllogfile = self._getlatestcpllog()
        cmd = os.path.join(self._case.get_value("SCRIPTSROOT"),"Tools","check_memory.pl")
        rc, out, err = run_cmd("%s -file1 %s -m 1.5"%(cmd, newestcpllogfile),ok_to_fail=True)
        if rc == 0:
            with open("TestStatus", "a") as fd:
                fd.write("PASS %s memleak\n"%(self._case.get_value("CASEBASEID")))
        else:
            with open(os.path.join(test_dir, "TestStatus.log"), "a") as fd:
                fd.write("memleak out: %s\n\nerror: %s"%(out,err))
            with open(os.path.join(test_dir, "TestStatus"), "a") as fd:
                fd.write("FAIL %s memleak\n"%(self._case.get_value("CASEBASEID")))

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

    def _getlatestcpllog(self):
        """
        find and return the latest cpl log file in the run directory
        """
        newestcpllogfile = min(glob.iglob(os.path.join(
                    self._case.get_value('RUNDIR'),'cpl.log.*')), key=os.path.getctime)
        return newestcpllogfile




class FakeTest(SystemTestsCommon):

    def fake_build(self, script, sharedlib_only=False, model_only=False):
        if (not sharedlib_only):
            exeroot = self._case.get_value("EXEROOT")
            cime_model = self._case.get_value("MODEL")
            modelexe = os.path.join(exeroot, "%s.exe"%cime_model)

            with open(modelexe, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(script)

            os.chmod(modelexe, 0755)
            self._case.set_value("BUILD_COMPLETE", True)
            self._case.flush()

class TESTRUNPASS(FakeTest):

    def build(self, sharedlib_only=False, model_only=False):
        rundir = self._case.get_value("RUNDIR")
        script = \
"""
echo Insta pass
echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID
""" % rundir
        self.fake_build(script,
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
cp %s/utils/python/tests/cpl.hi1.nc.test %s/%s.cpl.hi.0.nc.base
""" % (rundir, cimeroot, rundir, case)
        self.fake_build(script,
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
        self.fake_build(script,
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
        self.fake_build(script,
                        sharedlib_only=sharedlib_only, model_only=model_only)

