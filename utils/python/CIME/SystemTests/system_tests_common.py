"""
Base class for CIME system tests
"""
import shutil, glob
from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.XML.env_run import EnvRun
from CIME.utils import run_cmd
import CIME.build as build

class SystemTestsCommon(object):
    def __init__(self, caseroot=os.getcwd(), case=None, expectedrunvars=None):
        """
        initialize a CIME system test object, if the file LockedFiles/env_run.orig.xml
        does not exist copy the current env_run.xml file.  If it does exist restore values
        changed in a previous run of the test.
        """
        self._caseroot = caseroot
        self._case = case
        if os.path.isfile(os.path.join(caseroot, "LockedFiles", "env_run.orig.xml")):
            self.compare_env_run(expectedrunvars)
        elif os.path.isfile(os.path.join(caseroot, "env_run.xml")):
            lockedfiles = os.path.join(caseroot, "Lockedfiles")
            try:
                os.stat(lockedfiles)
            except:
                os.mkdir(lockedfiles)
            shutil.copy("env_run.xml",
                        os.path.join(lockedfiles, "env_run.orig.xml"))

    def build(self):
        build.case_build(self._caseroot, True)

    def run(self):
        run_cmd("case.run")
        return

    def report(self):
        return

    def check_mem_leak(self):
        """ TODO: incomplete """
        rundir = self._case.get_value("RUNDIR")
        cpllogfile = min(glob.iglob(os.path.join(rundir, "cpl.log*")), key=os.path.getctime)

    def compare_env_run(self, expected=None):
        f1obj = EnvRun(caseroot, "env_run.xml")
        f2obj = EnvRun(caseroot, os.path.join("LockedFiles", "env_run.orig.xml"))
        diffs = fi1obj.compare_xml(f2obj)
        for key in diffs.keys():
            if key in expected:
                logging.warn("  Resetting %s for test"%key)
                f1obj.set_value(key, f2obj.get_value(key, resolved=False))
            else:
                print "Found difference in %s: case: %s original value %s" %\
                    (key, diffs[key][0], diffs[key][1])
                print " Use option --force to run the test with this"\
                    " value or --reset to reset to original"
                return False
        return True

class TESTRUNPASS(SystemTestsCommon):
    def build(self):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        modelexe = os.path.join(exeroot, "%s.exe"%cime_model)
        with open(modelexe, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("echo Insta pass\n")
            f.write("echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID\n"%rundir)
        f.closed
        os.chmod(modelexe, 0755)
        self._case.set_value("BUILD_COMPLETE", "TRUE")
        self._case.flush()

class TESTRUNDIFF(SystemTestsCommon):
    def build(self):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        modelexe = os.path.join(exeroot, "%s.exe"%cime_model)
        with open(modelexe, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("echo Insta pass\n")
            f.write("echo SUCCESSFUL TERMINATION > %s/cpl.log.$LID\n"%rundir)
            f.write("cp %s/utils/python/tests/cpl.hi1.nc.test %s/%s.cpl.hi.0.nc.base\n"%
                    (cimeroot, rundir, case))
        f.closed
        os.chmod(modelexe, 0755)
        self._case.set_value("BUILD_COMPLETE", "TRUE")
        self._case.flush()

class TESTRUNFAIL(SystemTestsCommon):
    def build(self):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")
        rundir = self._case.get_value("RUNDIR")
        cimeroot = self._case.get_value("CIMEROOT")
        case = self._case.get_value("CASE")
        modelexe = os.path.join(exeroot, "%s.exe"%cime_model)
        with open(modelexe, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("echo Insta fail\n")
            f.write("echo model failed > %s/cpl.log.$LID\n"%rundir)
            f.write("$LID\n")
            f.write("exit -1\n")
        f.closed
        os.chmod(modelexe, 0755)
        self._case.set_value("BUILD_COMPLETE", "TRUE")
        self._case.flush()
