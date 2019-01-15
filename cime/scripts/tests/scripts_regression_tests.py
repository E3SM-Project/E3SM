#!/usr/bin/env python

"""
Script containing CIME python regression test suite. This suite should be run
to confirm overall CIME correctness.
"""

import glob, os, re, shutil, signal, sys, tempfile, \
    threading, time, logging, unittest, getpass, \
    filecmp, time

from xml.etree.ElementTree import ParseError

LIB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","lib")
sys.path.append(LIB_DIR)
# Remove all pyc files to ensure we're testing the right things
import subprocess, argparse
subprocess.call('/bin/rm -f $(find . -name "*.pyc")', shell=True, cwd=LIB_DIR)
import six
from six import assertRaisesRegex

import collections

from CIME.utils import run_cmd, run_cmd_no_fail, get_lids, get_current_commit, safe_copy
import get_tests
import CIME.test_scheduler, CIME.wait_for_tests
from  CIME.test_scheduler import TestScheduler
from  CIME.XML.compilers import Compilers
from  CIME.XML.env_run import EnvRun
from  CIME.XML.machines import Machines
from  CIME.XML.files import Files
from  CIME.case import Case
from  CIME.code_checker import check_code, get_all_checkable_files
from  CIME.test_status import *

SCRIPT_DIR  = CIME.utils.get_scripts_root()
TOOLS_DIR   = os.path.join(SCRIPT_DIR,"Tools")
TEST_COMPILER = None
GLOBAL_TIMEOUT = None
TEST_MPILIB = None
MACHINE     = None
FAST_ONLY   = False
NO_BATCH    = False
NO_CMAKE    = False
TEST_ROOT   = None
NO_TEARDOWN = False

os.environ["CIME_GLOBAL_WALLTIME"] = "0:05:00"

# pragma pylint: disable=protected-access
###############################################################################
def run_cmd_assert_result(test_obj, cmd, from_dir=None, expected_stat=0, env=None, verbose=False):
###############################################################################
    from_dir = os.getcwd() if from_dir is None else from_dir
    stat, output, errput = run_cmd(cmd, from_dir=from_dir, env=env, verbose=verbose)
    if expected_stat == 0:
        expectation = "SHOULD HAVE WORKED, INSTEAD GOT STAT %s" % stat
    else:
        expectation = "EXPECTED STAT %s, INSTEAD GOT STAT %s" % (expected_stat, stat)
    msg = \
"""
COMMAND: %s
FROM_DIR: %s
%s
OUTPUT: %s
ERRPUT: %s
""" % (cmd, from_dir, expectation, output, errput)
    test_obj.assertEqual(stat, expected_stat, msg=msg)

    return output

###############################################################################
def assert_test_status(test_obj, test_name, test_status_obj, test_phase, expected_stat):
###############################################################################
    test_status = test_status_obj.get_status(test_phase)
    test_obj.assertEqual(test_status, expected_stat, msg="Problem with {}: for phase '{}': has status '{}', expected '{}'".format(test_name, test_phase, test_status, expected_stat))

###############################################################################
class A_RunUnitTests(unittest.TestCase):
###############################################################################

    def test_resolve_variable_name(self):
        files = Files()
        machinefile = files.get_value("MACHINES_SPEC_FILE")
        self.assertTrue(os.path.isfile(machinefile),
                        msg="Path did not resolve to existing file %s" % machinefile)

    def test_unittests(self):
        # Finds all files contained in CIME/tests or its subdirectories that
        # match the pattern 'test*.py', and runs the unit tests found there
        # (i.e., tests defined using python's unittest module).
        #
        # This is analogous to running:
        #     python -m unittest discover -s CIME/tests -t .
        # from cime/scripts/lib
        #
        # Yes, that means we have a bunch of unit tests run from this one unit
        # test.

        testsuite = unittest.defaultTestLoader.discover(
            start_dir = os.path.join(LIB_DIR,"CIME","tests"),
            pattern = 'test*.py',
            top_level_dir = LIB_DIR)

        testrunner = unittest.TextTestRunner(buffer=False)

        # Disable logging; otherwise log messages written by code under test
        # clutter the unit test output
        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            results = testrunner.run(testsuite)
        finally:
            logging.getLogger().setLevel(log_lvl)

        self.assertTrue(results.wasSuccessful())

    def test_CIME_doctests(self):
        # Find and run all the doctests in the CIME directory tree
        run_cmd_assert_result(self, "PYTHONPATH=%s:$PYTHONPATH python -m doctest *.py 2>&1" % LIB_DIR, from_dir=os.path.join(LIB_DIR,"CIME"))

    def test_CIMEXML_doctests(self):
        # Find and run all the doctests in the XML directory tree
        run_cmd_assert_result(self, "PYTHONPATH=%s:$PYTHONPATH python -m doctest *.py 2>&1" % LIB_DIR, from_dir=os.path.join(LIB_DIR,"CIME","XML"))

###############################################################################
def make_fake_teststatus(path, testname, status, phase):
###############################################################################
    expect(phase in CORE_PHASES, "Bad phase '%s'" % phase)
    with TestStatus(test_dir=path, test_name=testname) as ts:
        for core_phase in CORE_PHASES:
            if core_phase == phase:
                ts.set_status(core_phase, status, comments=("time=42" if phase == RUN_PHASE else ""))
                break
            else:
                ts.set_status(core_phase, TEST_PASS_STATUS, comments=("time=42" if phase == RUN_PHASE else ""))

###############################################################################
def parse_test_status(line):
###############################################################################
    regex = re.compile(r"Test '(\w+)' finished with status '(\w+)'")
    m = regex.match(line)
    return m.groups()

###############################################################################
def kill_subprocesses(name=None, sig=signal.SIGKILL, expected_num_killed=None, tester=None):
###############################################################################
    # Kill all subprocesses
    proc_ids = CIME.utils.find_proc_id(proc_name=name, children_only=True)
    if (expected_num_killed is not None):
        tester.assertEqual(len(proc_ids), expected_num_killed,
                           msg="Expected to find %d processes to kill, found %d" % (expected_num_killed, len(proc_ids)))
    for proc_id in proc_ids:
        try:
            os.kill(proc_id, sig)
        except OSError:
            pass

###############################################################################
def kill_python_subprocesses(sig=signal.SIGKILL, expected_num_killed=None, tester=None):
###############################################################################
    kill_subprocesses("[Pp]ython", sig, expected_num_killed, tester)

###########################################################################
def assert_dashboard_has_build(tester, build_name, expected_count=1):
###########################################################################
    # Do not test E3SM dashboard if model is CESM
    if CIME.utils.get_model() == "e3sm":
        time.sleep(10) # Give chance for cdash to update

        wget_file = tempfile.mktemp()

        run_cmd_no_fail("wget https://my.cdash.org/index.php?project=ACME_test --no-check-certificate -O %s" % wget_file)

        raw_text = open(wget_file, "r").read()
        os.remove(wget_file)

        num_found = raw_text.count(build_name)
        tester.assertEqual(num_found, expected_count,
                           msg="Dashboard did not have expected num occurances of build name '%s'. Expected %s, found %s" % (build_name, expected_count, num_found))

###############################################################################
def setup_proxy():
###############################################################################
    if ("http_proxy" not in os.environ):
        proxy = MACHINE.get_value("PROXY")
        if (proxy is not None):
            os.environ["http_proxy"] = proxy
            return True

    return False

###############################################################################
class N_TestUnitTest(unittest.TestCase):
###############################################################################
    @classmethod
    def setUpClass(cls):
        cls._do_teardown = []
        cls._testroot = os.path.join(TEST_ROOT, 'TestUnitTests')
        cls._testdirs = []

    def _has_unit_test_support(self):
        if TEST_COMPILER is None:
            default_compiler = MACHINE.get_default_compiler()
            compiler = Compilers(MACHINE, compiler=default_compiler)
        else:
            compiler = Compilers(MACHINE, compiler=TEST_COMPILER)
        attrs = {'MPILIB': 'mpi-serial', 'compile_threaded': 'false'}
        pfunit_path = compiler.get_optional_compiler_node("PFUNIT_PATH",
                                                          attributes=attrs)
        if pfunit_path is None:
            return False
        else:
            return True

    def test_a_unit_test(self):
        cls = self.__class__
        if not self._has_unit_test_support():
            self.skipTest("Skipping TestUnitTest - PFUNIT_PATH not found for the default compiler on this machine")
        test_dir = os.path.join(cls._testroot,"unit_tester_test")
        cls._testdirs.append(test_dir)
        os.makedirs(test_dir)
        unit_test_tool = os.path.abspath(os.path.join(CIME.utils.get_cime_root(),"scripts","fortran_unit_testing","run_tests.py"))
        test_spec_dir = os.path.join(os.path.dirname(unit_test_tool),"Examples", "interpolate_1d", "tests")
        args = "--build-dir {} --test-spec-dir {}".format(test_dir, test_spec_dir)
        args += " --machine {}".format(MACHINE.get_machine_name())
        run_cmd_no_fail("{} {}".format(unit_test_tool, args))
        cls._do_teardown.append(test_dir)

    def test_b_cime_f90_unit_tests(self):
        cls = self.__class__
        if (FAST_ONLY):
            self.skipTest("Skipping slow test")

        if not self._has_unit_test_support():
            self.skipTest("Skipping TestUnitTest - PFUNIT_PATH not found for the default compiler on this machine")

        test_dir = os.path.join(cls._testroot,"driver_f90_tests")
        cls._testdirs.append(test_dir)
        os.makedirs(test_dir)
        test_spec_dir = CIME.utils.get_cime_root()
        unit_test_tool = os.path.abspath(os.path.join(test_spec_dir,"scripts","fortran_unit_testing","run_tests.py"))
        args = "--build-dir {} --test-spec-dir {}".format(test_dir, test_spec_dir)
        args += " --machine {}".format(MACHINE.get_machine_name())
        run_cmd_no_fail("{} {}".format(unit_test_tool, args))
        cls._do_teardown.append(test_dir)

    @classmethod
    def tearDownClass(cls):
        do_teardown = len(cls._do_teardown) > 0 and sys.exc_info() == (None, None, None) and not NO_TEARDOWN

        teardown_root = True
        for tfile in cls._testdirs:
            if tfile not in cls._do_teardown:
                print("Detected failed test or user request no teardown")
                print("Leaving case directory : %s"%tfile)
                teardown_root = False
            elif do_teardown:
                shutil.rmtree(tfile)

        if teardown_root and do_teardown:
            shutil.rmtree(cls._testroot)

###############################################################################
class J_TestCreateNewcase(unittest.TestCase):
###############################################################################
    @classmethod
    def setUpClass(cls):
        cls._testdirs = []
        cls._do_teardown = []
        cls._testroot = os.path.join(TEST_ROOT, 'TestCreateNewcase')

    def test_a_createnewcase(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'testcreatenewcase')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        args =  " --case %s --compset X --res f19_g16 --output-root %s --handle-preexisting-dirs=r" % (testdir, cls._testroot)
        if TEST_COMPILER is not None:
            args = args +  " --compiler %s"%TEST_COMPILER
        if TEST_MPILIB is not None:
            args = args +  " --mpilib %s"%TEST_MPILIB

        cls._testdirs.append(testdir)
        run_cmd_assert_result(self, "./create_newcase %s"%(args), from_dir=SCRIPT_DIR)
        self.assertTrue(os.path.exists(testdir))
        self.assertTrue(os.path.exists(os.path.join(testdir, "case.setup")))

        run_cmd_assert_result(self, "./case.setup", from_dir=testdir)
        run_cmd_assert_result(self, "./case.build", from_dir=testdir)

        with Case(testdir, read_only=False) as case:
            ntasks = case.get_value("NTASKS_ATM")
            case.set_value("NTASKS_ATM", ntasks+1)
        # this should fail with a locked file issue
        run_cmd_assert_result(self, "./case.build",
                              from_dir=testdir, expected_stat=1)

        run_cmd_assert_result(self, "./case.setup --reset", from_dir=testdir)
        run_cmd_assert_result(self, "./case.build", from_dir=testdir)
        with Case(testdir, read_only=False) as case:
            case.set_value("CHARGE_ACCOUNT", "fred")
        # this should not fail with a locked file issue
        run_cmd_assert_result(self, "./case.build",from_dir=testdir, expected_stat=0)

        run_cmd_assert_result(self, "./case.st_archive --test-all", from_dir=testdir)

        cls._do_teardown.append(testdir)

    def test_aa_no_flush_on_instantiate(self):
        testdir = os.path.join(self.__class__._testroot, 'testcreatenewcase')
        with Case(testdir, read_only=False) as case:
            for env_file in case._files:
                self.assertFalse(env_file.needsrewrite, msg="Instantiating a case should not trigger a flush call")

        with Case(testdir, read_only=False) as case:
            case.set_value("HIST_OPTION","nyears")
            runfile = case.get_env('run')
            self.assertTrue(runfile.needsrewrite, msg="Expected flush call not triggered")
            for env_file in case._files:
                if env_file != runfile:
                    self.assertFalse(env_file.needsrewrite, msg="Unexpected flush triggered for file {}"
                                     .format(env_file.filename))
            # Flush the file
            runfile.write()
            # set it again to the same value
            case.set_value("HIST_OPTION","nyears")
            # now the file should not need to be flushed
            for env_file in case._files:
                self.assertFalse(env_file.needsrewrite, msg="Unexpected flush triggered for file {}"
                                 .format(env_file.filename))

        # Check once more with a new instance
        with Case(testdir, read_only=False) as case:
            case.set_value("HIST_OPTION","nyears")
            for env_file in case._files:
                self.assertFalse(env_file.needsrewrite, msg="Unexpected flush triggered for file {}"
                                 .format(env_file.filename))

    def test_b_user_mods(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'testusermods')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)

        cls._testdirs.append(testdir)

        user_mods_dir = os.path.join(CIME.utils.get_python_libs_root(), "..", "tests", "user_mods_test1")
        args = " --case %s --compset X --res f19_g16 --user-mods-dir %s --output-root %s --handle-preexisting-dirs=r"% (testdir, user_mods_dir, cls._testroot)
        if TEST_COMPILER is not None:
            args = args + " --compiler %s"%TEST_COMPILER
        if TEST_MPILIB is not None:
            args = args +  " --mpilib %s"%TEST_MPILIB

        run_cmd_assert_result(self, "%s/create_newcase %s "
                              % (SCRIPT_DIR, args),from_dir=SCRIPT_DIR)

        self.assertTrue(os.path.isfile(os.path.join(testdir,"SourceMods","src.drv","somefile.F90")), msg="User_mods SourceMod missing")
        with open(os.path.join(testdir,"user_nl_cpl"),"r") as fd:
            contents = fd.read()
            self.assertTrue("a different cpl test option" in contents, msg="User_mods contents of user_nl_cpl missing")
            self.assertTrue("a cpl namelist option" in contents, msg="User_mods contents of user_nl_cpl missing")
        cls._do_teardown.append(testdir)

    def test_c_create_clone_keepexe(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'test_create_clone_keepexe')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        prevtestdir = cls._testdirs[0]
        user_mods_dir = os.path.join(CIME.utils.get_python_libs_root(), "..", "tests", "user_mods_test3")

        cmd = "%s/create_clone --clone %s --case %s --keepexe --user-mods-dir %s" \
              % (SCRIPT_DIR, prevtestdir, testdir, user_mods_dir)
        run_cmd_assert_result(self, cmd, from_dir=SCRIPT_DIR, expected_stat=1)

    def test_d_create_clone_new_user(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'test_create_clone_new_user')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        prevtestdir = cls._testdirs[0]
        cls._testdirs.append(testdir)
        # change the USER and CIME_OUTPUT_ROOT to nonsense values
        # this is intended as a test of whether create_clone is independent of user
        run_cmd_assert_result(self, "./xmlchange USER=this_is_not_a_user",
                              from_dir=prevtestdir)

        fakeoutputroot = cls._testroot.replace(os.environ.get("USER"), "this_is_not_a_user")
        run_cmd_assert_result(self, "./xmlchange CIME_OUTPUT_ROOT=%s"%fakeoutputroot,
                              from_dir=prevtestdir)

        # this test should pass (user name is replaced)
        run_cmd_assert_result(self, "%s/create_clone --clone %s --case %s " %
                              (SCRIPT_DIR, prevtestdir, testdir),from_dir=SCRIPT_DIR)

        shutil.rmtree(testdir)
        # this test should pass
        run_cmd_assert_result(self, "%s/create_clone --clone %s --case %s --cime-output-root %s" %
                              (SCRIPT_DIR, prevtestdir, testdir, cls._testroot),from_dir=SCRIPT_DIR)

        cls._do_teardown.append(testdir)

    def test_e_xmlquery(self):
        # Set script and script path
        xmlquery = "./xmlquery"
        cls = self.__class__
        casedir  =  cls._testdirs[0]

        # Check for environment
        self.assertTrue(os.path.isdir(SCRIPT_DIR))
        self.assertTrue(os.path.isdir(TOOLS_DIR))
        self.assertTrue(os.path.isfile(os.path.join(casedir,xmlquery)))

        # Test command line options
        with Case(casedir, read_only=True) as case:
            STOP_N = case.get_value("STOP_N")
            COMP_CLASSES = case.get_values("COMP_CLASSES")
            BUILD_COMPLETE = case.get_value("BUILD_COMPLETE")
            cmd = xmlquery + " STOP_N --value"
            output = run_cmd_no_fail(cmd, from_dir=casedir)
            self.assertTrue(output == str(STOP_N), msg="%s != %s"%(output, STOP_N))
            cmd = xmlquery + " BUILD_COMPLETE --value"
            output = run_cmd_no_fail(cmd, from_dir=casedir)
            self.assertTrue(output == "TRUE", msg="%s != %s"%(output, BUILD_COMPLETE))
            # we expect DOCN_MODE to be undefined in this X compset
            # this test assures that we do not try to resolve this as a compvar
            cmd = xmlquery + " DOCN_MODE --value"
            _, output, error = run_cmd(cmd, from_dir=casedir)
            self.assertTrue(error == "ERROR:  No results found for variable DOCN_MODE",
                            msg="unexpected result for DOCN_MODE, output {}, error {}".
                            format(output, error))

            for comp in COMP_CLASSES:
                caseresult = case.get_value("NTASKS_%s"%comp)
                cmd = xmlquery + " NTASKS_%s --value"%comp
                output = run_cmd_no_fail(cmd, from_dir=casedir)
                self.assertTrue(output == str(caseresult), msg="%s != %s"%(output, caseresult))
                cmd = xmlquery + " NTASKS --subgroup %s --value"%comp
                output = run_cmd_no_fail(cmd, from_dir=casedir)
                self.assertTrue(output == str(caseresult), msg="%s != %s"%(output, caseresult))
            if MACHINE.has_batch_system():
                JOB_QUEUE = case.get_value("JOB_QUEUE", subgroup="case.run")
                cmd = xmlquery + " JOB_QUEUE --subgroup case.run --value"
                output = run_cmd_no_fail(cmd, from_dir=casedir)
                self.assertTrue(output == JOB_QUEUE, msg="%s != %s"%(output, JOB_QUEUE))

            cmd = xmlquery + " --listall"
            run_cmd_no_fail(cmd, from_dir=casedir)

        cls._do_teardown.append(cls._testroot)

    def test_f_createnewcase_with_user_compset(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'testcreatenewcase_with_user_compset')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)

        cls._testdirs.append(testdir)

        pesfile = os.path.join("..","src","drivers","mct","cime_config","config_pes.xml")
        args =  "--case %s --compset 2000_SATM_XLND_SICE_SOCN_XROF_XGLC_SWAV  --pesfile %s --res f19_g16 --output-root %s --handle-preexisting-dirs=r" % (testdir, pesfile, cls._testroot)
        if CIME.utils.get_model() == "cesm":
            args += " --run-unsupported"
        if TEST_COMPILER is not None:
            args += " --compiler %s"%TEST_COMPILER
        if TEST_MPILIB is not None:
            args = args +  " --mpilib %s"%TEST_MPILIB

        run_cmd_assert_result(self, "%s/create_newcase %s"%(SCRIPT_DIR, args), from_dir=SCRIPT_DIR)
        run_cmd_assert_result(self, "./case.setup", from_dir=testdir)
        run_cmd_assert_result(self, "./case.build", from_dir=testdir)

        cls._do_teardown.append(testdir)

    def test_g_createnewcase_with_user_compset_and_env_mach_pes(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'testcreatenewcase_with_user_compset_and_env_mach_pes')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        previous_testdir = cls._testdirs[-1]
        cls._testdirs.append(testdir)

        pesfile = os.path.join(previous_testdir,"env_mach_pes.xml")
        args =  "--case %s --compset 2000_SATM_XLND_SICE_SOCN_XROF_XGLC_SWAV --pesfile %s --res f19_g16 --output-root %s --handle-preexisting-dirs=r" % (testdir, pesfile, cls._testroot)
        if CIME.utils.get_model() == "cesm":
            args += " --run-unsupported"
        if TEST_COMPILER is not None:
            args += " --compiler %s"%TEST_COMPILER
        if TEST_MPILIB is not None:
            args +=  " --mpilib %s"%TEST_MPILIB

        run_cmd_assert_result(self, "%s/create_newcase %s"%(SCRIPT_DIR, args), from_dir=SCRIPT_DIR)
        run_cmd_assert_result(self, "diff env_mach_pes.xml %s"%(previous_testdir), from_dir=testdir)
        # this line should cause the diff to fail (I assume no machine is going to default to 17 tasks)
        run_cmd_assert_result(self, "./xmlchange NTASKS=17", from_dir=testdir)
        run_cmd_assert_result(self, "diff env_mach_pes.xml %s"%(previous_testdir), from_dir=testdir,
                              expected_stat=1)

        cls._do_teardown.append(testdir)

    def test_h_primary_component(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'testprimarycomponent')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)

        cls._testdirs.append(testdir)
        args = " --case CreateNewcaseTest --script-root %s --compset X --res f19_g16 --output-root %s --handle-preexisting-dirs u" % (testdir, cls._testroot)
        if TEST_COMPILER is not None:
            args += " --compiler %s"%TEST_COMPILER
        if TEST_MPILIB is not None:
            args +=  " --mpilib %s"%TEST_MPILIB

        run_cmd_assert_result(self, "%s/create_newcase %s" % (SCRIPT_DIR, args), from_dir=SCRIPT_DIR)
        self.assertTrue(os.path.exists(testdir))
        self.assertTrue(os.path.exists(os.path.join(testdir, "case.setup")))

        with Case(testdir, read_only=False) as case:
            case._compsetname = case.get_value("COMPSET")
            case.set_comp_classes(case.get_values("COMP_CLASSES"))
            primary = case._find_primary_component()
            self.assertEqual(primary, "drv", msg="primary component test expected drv but got %s"%primary)
            # now we are going to corrupt the case so that we can do more primary_component testing
            case.set_valid_values("COMP_GLC","%s,fred"%case.get_value("COMP_GLC"))
            case.set_value("COMP_GLC","fred")
            primary = case._find_primary_component()
            self.assertEqual(primary, "fred", msg="primary component test expected fred but got %s"%primary)
            case.set_valid_values("COMP_ICE","%s,wilma"%case.get_value("COMP_ICE"))
            case.set_value("COMP_ICE","wilma")
            primary = case._find_primary_component()
            self.assertEqual(primary, "wilma", msg="primary component test expected wilma but got %s"%primary)

            case.set_valid_values("COMP_OCN","%s,bambam,docn"%case.get_value("COMP_OCN"))
            case.set_value("COMP_OCN","bambam")
            primary = case._find_primary_component()
            self.assertEqual(primary, "bambam", msg="primary component test expected bambam but got %s"%primary)

            case.set_valid_values("COMP_LND","%s,barney"%case.get_value("COMP_LND"))
            case.set_value("COMP_LND","barney")
            primary = case._find_primary_component()
            # This is a "J" compset
            self.assertEqual(primary, "allactive", msg="primary component test expected allactive but got %s"%primary)
            case.set_value("COMP_OCN","docn")
            case.set_valid_values("COMP_LND","%s,barney"%case.get_value("COMP_LND"))
            case.set_value("COMP_LND","barney")
            primary = case._find_primary_component()
            self.assertEqual(primary, "barney", msg="primary component test expected barney but got %s"%primary)
            case.set_valid_values("COMP_ATM","%s,wilma"%case.get_value("COMP_ATM"))
            case.set_value("COMP_ATM","wilma")
            primary = case._find_primary_component()
            self.assertEqual(primary, "wilma", msg="primary component test expected wilma but got %s"%primary)
            # this is a "E" compset
            case._compsetname = case._compsetname.replace("XOCN","DOCN%SOM")
            primary = case._find_primary_component()
            self.assertEqual(primary, "allactive", msg="primary component test expected allactive but got %s"%primary)
            # finally a "B" compset
            case.set_value("COMP_OCN","bambam")
            primary = case._find_primary_component()
            self.assertEqual(primary, "allactive", msg="primary component test expected allactive but got %s"%primary)

        cls._do_teardown.append(testdir)

    def test_j_createnewcase_user_compset_vs_alias(self):
        """
        Create a compset using the alias and another compset using the full compset name
        and make sure they are the same by comparing the namelist files in CaseDocs.
        Ignore the modelio files and clean the directory names out first.
        """

        cls = self.__class__

        testdir1 = os.path.join(cls._testroot, 'testcreatenewcase_user_compset')
        if os.path.exists(testdir1):
            shutil.rmtree(testdir1)
        cls._testdirs.append(testdir1)
        args = ' --case CreateNewcaseTest --script-root {} --compset 2000_DATM%NYF_SLND_SICE_DOCN%SOMAQP_SROF_SGLC_SWAV --res f19_g16 --output-root {} --handle-preexisting-dirs u' .format(testdir1, cls._testroot)
        if CIME.utils.get_model() == "cesm":
            args += " --run-unsupported"
        if TEST_COMPILER is not None:
            args += " --compiler %s"%TEST_COMPILER
        if TEST_MPILIB is not None:
            args +=  " --mpilib %s"%TEST_MPILIB

        run_cmd_assert_result(self, "{}/create_newcase {}" .format (SCRIPT_DIR, args), from_dir=SCRIPT_DIR)
        run_cmd_assert_result(self, "./case.setup ", from_dir=testdir1)
        run_cmd_assert_result(self, "./preview_namelists ", from_dir=testdir1)

        dir1 = os.path.join(testdir1,"CaseDocs")
        dir2 = os.path.join(testdir1,"CleanCaseDocs")
        os.mkdir(dir2)
        for _file in os.listdir(dir1):
            if "modelio" in _file:
                continue
            with open(os.path.join(dir1,_file),"r") as fi:
                file_text = fi.read()
                file_text = file_text.replace(os.path.basename(testdir1),"PATH")
            with open(os.path.join(dir2,_file), "w") as fo:
                fo.write(file_text)
        cleancasedocs1 = dir2

        testdir2 = os.path.join(cls._testroot, 'testcreatenewcase_alias_compset')
        if os.path.exists(testdir2):
            shutil.rmtree(testdir2)
        cls._testdirs.append(testdir2)
        args = ' --case CreateNewcaseTest --script-root {} --compset ADSOMAQP --res f19_g16 --output-root {} --handle-preexisting-dirs u'.format(testdir2, cls._testroot)
        if CIME.utils.get_model() == "cesm":
            args += " --run-unsupported"
        if TEST_COMPILER is not None:
            args += " --compiler %s"%TEST_COMPILER
        if TEST_MPILIB is not None:
            args +=  " --mpilib %s"%TEST_MPILIB

        run_cmd_assert_result(self, "{}/create_newcase {}".format(SCRIPT_DIR, args), from_dir=SCRIPT_DIR)
        run_cmd_assert_result(self, "./case.setup ", from_dir=testdir2)
        run_cmd_assert_result(self, "./preview_namelists ", from_dir=testdir2)

        dir1 = os.path.join(testdir2,"CaseDocs")
        dir2 = os.path.join(testdir2,"CleanCaseDocs")
        os.mkdir(dir2)
        for _file in os.listdir(dir1):
            if "modelio" in _file:
                continue
            with open(os.path.join(dir1,_file),"r") as fi:
                file_text = fi.read()
                file_text = file_text.replace(os.path.basename(testdir2),"PATH")
            with open(os.path.join(dir2,_file), "w") as fo:
                fo.write(file_text)

        cleancasedocs2 = dir2
        dcmp = filecmp.dircmp(cleancasedocs1, cleancasedocs2)
        self.assertTrue(len(dcmp.diff_files) == 0, "CaseDocs differ {}".format(dcmp.diff_files))

        cls._do_teardown.append(testdir1)
        cls._do_teardown.append(testdir2)


    def test_k_append_config(self):
        machlist_before = MACHINE.list_available_machines()
        self.assertEqual(len(machlist_before)>1, True, msg="Problem reading machine list")

        newmachfile = os.path.join(CIME.utils.get_cime_root(),"config",
                                   "xml_schemas","config_machines_template.xml")
        MACHINE.read(newmachfile)
        machlist_after = MACHINE.list_available_machines()

        self.assertEqual(len(machlist_after)-len(machlist_before), 1, msg="Not able to append config_machines.xml {} {}".format(len(machlist_after), len(machlist_before)))
        self.assertEqual("mymachine" in machlist_after, True, msg="Not able to append config_machines.xml")


    def test_m_createnewcase_alternate_drivers(self):
        # Test that case.setup runs for nuopc and moab drivers
        cls = self.__class__
        for driver in ("nuopc", "moab"):
            testdir = os.path.join(cls._testroot, 'testcreatenewcase.{}'.format( driver))
            if os.path.exists(testdir):
                shutil.rmtree(testdir)
            args =  " --driver {} --case {} --compset X --res f19_g16 --output-root {} --handle-preexisting-dirs=r".format(driver, testdir, cls._testroot)
            if CIME.utils.get_model() == "cesm":
                args += " --run-unsupported"
            if TEST_COMPILER is not None:
                args = args +  " --compiler %s"%TEST_COMPILER
            if TEST_MPILIB is not None:
                args = args +  " --mpilib %s"%TEST_MPILIB

            cls._testdirs.append(testdir)
            run_cmd_assert_result(self, "./create_newcase %s"%(args), from_dir=SCRIPT_DIR)
            self.assertTrue(os.path.exists(testdir))
            self.assertTrue(os.path.exists(os.path.join(testdir, "case.setup")))

            run_cmd_assert_result(self, "./case.setup", from_dir=testdir)
            with Case(testdir, read_only=False) as case:
                comp_interface = case.get_value("COMP_INTERFACE")
                self.assertTrue(driver == comp_interface, msg="%s != %s"%(driver, comp_interface))

            cls._do_teardown.append(testdir)


    @classmethod
    def tearDownClass(cls):
        do_teardown = len(cls._do_teardown) > 0 and sys.exc_info() == (None, None, None) and not NO_TEARDOWN

        for tfile in cls._testdirs:
            if tfile not in cls._do_teardown:
                print("Detected failed test or user request no teardown")
                print("Leaving case directory : %s"%tfile)
            elif do_teardown:
                try:
                    print ("Attempt to remove directory {}".format(tfile))
                    shutil.rmtree(tfile)
                except BaseException:
                    print("Could not remove directory {}".format(tfile))

###############################################################################
class M_TestWaitForTests(unittest.TestCase):
###############################################################################

    ###########################################################################
    def setUp(self):
    ###########################################################################
        self._testroot = os.path.join(TEST_ROOT,"TestWaitForTests")

        # basic tests
        self._testdir_all_pass    = os.path.join(self._testroot, 'scripts_regression_tests.testdir_all_pass')
        self._testdir_with_fail   = os.path.join(self._testroot, 'scripts_regression_tests.testdir_with_fail')
        self._testdir_unfinished  = os.path.join(self._testroot, 'scripts_regression_tests.testdir_unfinished')
        self._testdir_unfinished2 = os.path.join(self._testroot, 'scripts_regression_tests.testdir_unfinished2')

        # live tests
        self._testdir_teststatus1 = os.path.join(self._testroot, 'scripts_regression_tests.testdir_teststatus1')
        self._testdir_teststatus2 = os.path.join(self._testroot, 'scripts_regression_tests.testdir_teststatus2')

        self._testdirs = [self._testdir_all_pass, self._testdir_with_fail, self._testdir_unfinished, self._testdir_unfinished2,
                          self._testdir_teststatus1, self._testdir_teststatus2]
        basic_tests = self._testdirs[:self._testdirs.index(self._testdir_teststatus1)]

        for testdir in self._testdirs:
            if os.path.exists(testdir):
                shutil.rmtree(testdir)
            os.makedirs(testdir)

        for r in range(10):
            for testdir in basic_tests:
                os.makedirs(os.path.join(testdir, str(r)))
                make_fake_teststatus(os.path.join(testdir, str(r)), "Test_%d" % r, TEST_PASS_STATUS, RUN_PHASE)

        make_fake_teststatus(os.path.join(self._testdir_with_fail,   "5"), "Test_5", TEST_FAIL_STATUS, RUN_PHASE)
        make_fake_teststatus(os.path.join(self._testdir_unfinished,  "5"), "Test_5", TEST_PEND_STATUS, RUN_PHASE)
        make_fake_teststatus(os.path.join(self._testdir_unfinished2, "5"), "Test_5", TEST_PASS_STATUS, SUBMIT_PHASE)

        integration_tests = self._testdirs[len(basic_tests):]
        for integration_test in integration_tests:
            os.makedirs(os.path.join(integration_test, "0"))
            make_fake_teststatus(os.path.join(integration_test, "0"), "Test_0", TEST_PASS_STATUS, CORE_PHASES[0])

        # Set up proxy if possible
        self._unset_proxy = setup_proxy()

        self._thread_error = None

    ###########################################################################
    def tearDown(self):
    ###########################################################################
        do_teardown = sys.exc_info() == (None, None, None) and not NO_TEARDOWN

        if do_teardown:
            for testdir in self._testdirs:
                shutil.rmtree(testdir)

        kill_subprocesses()

        if (self._unset_proxy):
            del os.environ["http_proxy"]

    ###########################################################################
    def simple_test(self, testdir, expected_results, extra_args="", build_name=None):
    ###########################################################################
        # Need these flags to test dashboard if e3sm
        if CIME.utils.get_model() == "e3sm" and build_name is not None:
            extra_args += " -b %s" % build_name

        expected_stat = 0 if expected_results == ["PASS"]*len(expected_results) else CIME.utils.TESTS_FAILED_ERR_CODE
        output = run_cmd_assert_result(self, "%s/wait_for_tests -p ACME_test */TestStatus %s" % (TOOLS_DIR, extra_args),
                                       from_dir=testdir, expected_stat=expected_stat)

        lines = [line for line in output.splitlines() if line.startswith("Test '")]
        self.assertEqual(len(lines), len(expected_results))
        for idx, line in enumerate(lines):
            testname, status = parse_test_status(line)
            self.assertEqual(status, expected_results[idx])
            self.assertEqual(testname, "Test_%d" % idx)

    ###########################################################################
    def threaded_test(self, testdir, expected_results, extra_args="", build_name=None):
    ###########################################################################
        try:
            self.simple_test(testdir, expected_results, extra_args, build_name)
        except AssertionError as e:
            self._thread_error = str(e)

    ###########################################################################
    def test_wait_for_test_all_pass(self):
    ###########################################################################
        self.simple_test(self._testdir_all_pass, ["PASS"] * 10)

    ###########################################################################
    def test_wait_for_test_with_fail(self):
    ###########################################################################
        expected_results = ["FAIL" if item == 5 else "PASS" for item in range(10)]
        self.simple_test(self._testdir_with_fail, expected_results)

    ###########################################################################
    def test_wait_for_test_no_wait(self):
    ###########################################################################
        expected_results = ["PEND" if item == 5 else "PASS" for item in range(10)]
        self.simple_test(self._testdir_unfinished, expected_results, "-n")

    ###########################################################################
    def test_wait_for_test_timeout(self):
    ###########################################################################
        expected_results = ["PEND" if item == 5 else "PASS" for item in range(10)]
        self.simple_test(self._testdir_unfinished, expected_results, "--timeout=3")

    ###########################################################################
    def test_wait_for_test_wait_for_pend(self):
    ###########################################################################
        run_thread = threading.Thread(target=self.threaded_test, args=(self._testdir_unfinished, ["PASS"] * 10))
        run_thread.daemon = True
        run_thread.start()

        time.sleep(5) # Kinda hacky

        self.assertTrue(run_thread.isAlive(), msg="wait_for_tests should have waited")

        with TestStatus(test_dir=os.path.join(self._testdir_unfinished, "5")) as ts:
            ts.set_status(RUN_PHASE, TEST_PASS_STATUS)

        run_thread.join(timeout=10)

        self.assertFalse(run_thread.isAlive(), msg="wait_for_tests should have finished")

        self.assertTrue(self._thread_error is None, msg="Thread had failure: %s" % self._thread_error)

    ###########################################################################
    def test_wait_for_test_wait_for_missing_run_phase(self):
    ###########################################################################
        run_thread = threading.Thread(target=self.threaded_test, args=(self._testdir_unfinished2, ["PASS"] * 10))
        run_thread.daemon = True
        run_thread.start()

        time.sleep(5) # Kinda hacky

        self.assertTrue(run_thread.isAlive(), msg="wait_for_tests should have waited")

        with TestStatus(test_dir=os.path.join(self._testdir_unfinished2, "5")) as ts:
            ts.set_status(RUN_PHASE, TEST_PASS_STATUS)

        run_thread.join(timeout=10)

        self.assertFalse(run_thread.isAlive(), msg="wait_for_tests should have finished")

        self.assertTrue(self._thread_error is None, msg="Thread had failure: %s" % self._thread_error)

    ###########################################################################
    def test_wait_for_test_wait_kill(self):
    ###########################################################################
        expected_results = ["PEND" if item == 5 else "PASS" for item in range(10)]
        run_thread = threading.Thread(target=self.threaded_test, args=(self._testdir_unfinished, expected_results))
        run_thread.daemon = True
        run_thread.start()

        time.sleep(5)

        self.assertTrue(run_thread.isAlive(), msg="wait_for_tests should have waited")

        kill_python_subprocesses(signal.SIGTERM, expected_num_killed=1, tester=self)

        run_thread.join(timeout=10)

        self.assertFalse(run_thread.isAlive(), msg="wait_for_tests should have finished")

        self.assertTrue(self._thread_error is None, msg="Thread had failure: %s" % self._thread_error)

    ###########################################################################
    def test_wait_for_test_cdash_pass(self):
    ###########################################################################
        expected_results = ["PASS"] * 10
        run_thread = threading.Thread(target=self.threaded_test,
                                      args=(self._testdir_all_pass, expected_results, "", "regression_test_pass"))
        run_thread.daemon = True
        run_thread.start()

        run_thread.join(timeout=10)

        self.assertFalse(run_thread.isAlive(), msg="wait_for_tests should have finished")

        self.assertTrue(self._thread_error is None, msg="Thread had failure: %s" % self._thread_error)

        assert_dashboard_has_build(self, "regression_test_pass")

    ###########################################################################
    def test_wait_for_test_cdash_kill(self):
    ###########################################################################
        expected_results = ["PEND" if item == 5 else "PASS" for item in range(10)]
        run_thread = threading.Thread(target=self.threaded_test,
                                      args=(self._testdir_unfinished, expected_results, "", "regression_test_kill"))
        run_thread.daemon = True
        run_thread.start()

        time.sleep(5)

        self.assertTrue(run_thread.isAlive(), msg="wait_for_tests should have waited")

        kill_python_subprocesses(signal.SIGTERM, expected_num_killed=1, tester=self)

        run_thread.join(timeout=10)

        self.assertFalse(run_thread.isAlive(), msg="wait_for_tests should have finished")
        self.assertTrue(self._thread_error is None, msg="Thread had failure: %s" % self._thread_error)

        assert_dashboard_has_build(self, "regression_test_kill")

        if CIME.utils.get_model() == "e3sm":
            cdash_result_dir = os.path.join(self._testdir_unfinished, "Testing")
            tag_file         = os.path.join(cdash_result_dir, "TAG")
            self.assertTrue(os.path.isdir(cdash_result_dir))
            self.assertTrue(os.path.isfile(tag_file))

            tag = open(tag_file, "r").readlines()[0].strip()
            xml_file = os.path.join(cdash_result_dir, tag, "Test.xml")
            self.assertTrue(os.path.isfile(xml_file))

            xml_contents = open(xml_file, "r").read()
            self.assertTrue(r'<TestList><Test>Test_0</Test><Test>Test_1</Test><Test>Test_2</Test><Test>Test_3</Test><Test>Test_4</Test><Test>Test_5</Test><Test>Test_6</Test><Test>Test_7</Test><Test>Test_8</Test><Test>Test_9</Test></TestList>'
                            in xml_contents)
            self.assertTrue(r'<Test Status="notrun"><Name>Test_5</Name>' in xml_contents)

            # TODO: Any further checking of xml output worth doing?

    ###########################################################################
    def live_test_impl(self, testdir, expected_results, last_phase, last_status):
    ###########################################################################
        run_thread = threading.Thread(target=self.threaded_test, args=(testdir, expected_results))
        run_thread.daemon = True
        run_thread.start()

        time.sleep(5)

        self.assertTrue(run_thread.isAlive(), msg="wait_for_tests should have waited")

        for core_phase in CORE_PHASES[1:]:
            with TestStatus(test_dir=os.path.join(self._testdir_teststatus1, "0")) as ts:
                ts.set_status(core_phase, last_status if core_phase == last_phase else TEST_PASS_STATUS)

            time.sleep(5)

            if core_phase != last_phase:
                self.assertTrue(run_thread.isAlive(), msg="wait_for_tests should have waited after passing phase {}".format(core_phase))
            else:
                run_thread.join(timeout=10)
                self.assertFalse(run_thread.isAlive(), msg="wait_for_tests should have finished after phase {}".format(core_phase))
                break

        self.assertTrue(self._thread_error is None, msg="Thread had failure: %s" % self._thread_error)

    ###########################################################################
    def test_wait_for_test_test_status_integration_pass(self):
    ###########################################################################
        self.live_test_impl(self._testdir_teststatus1, ["PASS"], RUN_PHASE, TEST_PASS_STATUS)

    ###########################################################################
    def test_wait_for_test_test_status_integration_submit_fail(self):
    ###########################################################################
        self.live_test_impl(self._testdir_teststatus1, ["FAIL"], SUBMIT_PHASE, TEST_FAIL_STATUS)

###############################################################################
class TestCreateTestCommon(unittest.TestCase):
###############################################################################

    ###########################################################################
    def setUp(self):
    ###########################################################################
        self._thread_error      = None
        self._unset_proxy       = setup_proxy()
        self._machine           = MACHINE.get_machine_name()
        self._compiler          = MACHINE.get_default_compiler() if TEST_COMPILER is None else TEST_COMPILER
        self._baseline_name     = "fake_testing_only_%s" % CIME.utils.get_timestamp()
        self._baseline_area     = os.path.join(TEST_ROOT, "baselines")
        self._testroot          = TEST_ROOT
        self._hasbatch          = MACHINE.has_batch_system() and not NO_BATCH
        self._do_teardown       = not NO_TEARDOWN

    ###########################################################################
    def tearDown(self):
    ###########################################################################
        kill_subprocesses()

        if (self._unset_proxy):
            del os.environ["http_proxy"]

        files_to_clean = []

        baselines = os.path.join(self._baseline_area, self._baseline_name)
        if (os.path.isdir(baselines)):
            files_to_clean.append(baselines)

        for test_id in ["master", self._baseline_name]:
            for leftover in glob.glob(os.path.join(self._testroot, "*%s*" % test_id)):
                files_to_clean.append(leftover)

        do_teardown = self._do_teardown and sys.exc_info() == (None, None, None)
        if (not do_teardown):
            print("Detected failed test or user request no teardown")
            print("Leaving files:")
            for file_to_clean in files_to_clean:
                print(" " + file_to_clean)
        else:
            # For batch machines need to avoid race condition as batch system
            # finishes I/O for the case.
            if self._hasbatch:
                time.sleep(5)

            for file_to_clean in files_to_clean:
                if (os.path.isdir(file_to_clean)):
                    shutil.rmtree(file_to_clean)
                else:
                    os.remove(file_to_clean)

    ###########################################################################
    def _create_test(self, extra_args, test_id=None, pre_run_errors=False, run_errors=False, env_changes=""):
    ###########################################################################
        test_id = CIME.utils.get_timestamp() if test_id is None else test_id
        extra_args.append("-t {}".format(test_id))
        extra_args.append("--baseline-root {}".format(self._baseline_area))
        if NO_BATCH:
            extra_args.append("--no-batch")
        if TEST_COMPILER and ([extra_arg for extra_arg in extra_args if "--compiler" in extra_arg] == []):
            extra_args.append("--compiler={}".format(TEST_COMPILER))
        if TEST_MPILIB and ([extra_arg for extra_arg in extra_args if "--mpilib" in extra_arg] == []):
            extra_args.append("--mpilib={}".format(TEST_MPILIB))
        extra_args.append("--test-root={0} --output-root={0}".format(TEST_ROOT))

        full_run = (set(extra_args) & set(["-n", "--namelist-only", "--no-setup", "--no-build"])) == set()

        if self._hasbatch:
            expected_stat = 0 if not pre_run_errors else CIME.utils.TESTS_FAILED_ERR_CODE
        else:
            expected_stat = 0 if not pre_run_errors and not run_errors else CIME.utils.TESTS_FAILED_ERR_CODE

        run_cmd_assert_result(self, "{} {}/create_test {}".format(env_changes, SCRIPT_DIR, " ".join(extra_args)),
                              expected_stat=expected_stat)

        if full_run:
            self._wait_for_tests(test_id, expect_works=(not pre_run_errors and not run_errors))

    ###########################################################################
    def _wait_for_tests(self, test_id, expect_works=True):
    ###########################################################################
        if self._hasbatch:
            timeout_arg = "--timeout={}".format(GLOBAL_TIMEOUT) if GLOBAL_TIMEOUT is not None else ""
            expected_stat = 0 if expect_works else CIME.utils.TESTS_FAILED_ERR_CODE
            run_cmd_assert_result(self, "{}/wait_for_tests {} *{}/TestStatus".format(TOOLS_DIR, timeout_arg, test_id),
                                  from_dir=self._testroot, expected_stat=expected_stat)

###############################################################################
class O_TestTestScheduler(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_a_phases(self):
    ###########################################################################
        # exclude the MEMLEAK tests here.
        tests = get_tests.get_full_test_names(["cime_test_only",
                                                       "^TESTMEMLEAKFAIL_P1.f09_g16.X",
                                                       "^TESTMEMLEAKPASS_P1.f09_g16.X",
                                                       "^TESTRUNSTARCFAIL_P1.f19_g16_rx1.A",
                                                       "^TESTTESTDIFF_P1.f19_g16_rx1.A",
                                                       "^TESTBUILDFAILEXC_P1.f19_g16_rx1.A",
                                                       "^TESTRUNFAILEXC_P1.f19_g16_rx1.A"],
                                                      self._machine, self._compiler)
        self.assertEqual(len(tests), 3)
        ct = TestScheduler(tests, test_root=TEST_ROOT, output_root=TEST_ROOT,
                           compiler=self._compiler, mpilib=TEST_MPILIB)

        build_fail_test = [item for item in tests if "TESTBUILDFAIL" in item][0]
        run_fail_test   = [item for item in tests if "TESTRUNFAIL" in item][0]
        pass_test       = [item for item in tests if "TESTRUNPASS" in item][0]

        self.assertTrue("BUILDFAIL" in build_fail_test, msg="Wrong test '%s'" % build_fail_test)
        self.assertTrue("RUNFAIL" in run_fail_test, msg="Wrong test '%s'" % run_fail_test)
        self.assertTrue("RUNPASS" in pass_test, msg="Wrong test '%s'" % pass_test)

        for idx, phase in enumerate(ct._phases):
            for test in ct._tests:
                if (phase == CIME.test_scheduler.TEST_START):
                    continue
                elif (phase == MODEL_BUILD_PHASE):
                    ct._update_test_status(test, phase, TEST_PEND_STATUS)

                    if (test == build_fail_test):
                        ct._update_test_status(test, phase, TEST_FAIL_STATUS)
                        self.assertTrue(ct._is_broken(test))
                        self.assertFalse(ct._work_remains(test))
                    else:
                        ct._update_test_status(test, phase, TEST_PASS_STATUS)
                        self.assertFalse(ct._is_broken(test))
                        self.assertTrue(ct._work_remains(test))

                elif (phase == RUN_PHASE):
                    if (test == build_fail_test):
                        with self.assertRaises(SystemExit):
                            ct._update_test_status(test, phase, TEST_PEND_STATUS)
                    else:
                        ct._update_test_status(test, phase, TEST_PEND_STATUS)
                        self.assertFalse(ct._work_remains(test))

                        if (test == run_fail_test):
                            ct._update_test_status(test, phase, TEST_FAIL_STATUS)
                            self.assertTrue(ct._is_broken(test))
                        else:
                            ct._update_test_status(test, phase, TEST_PASS_STATUS)
                            self.assertFalse(ct._is_broken(test))

                    self.assertFalse(ct._work_remains(test))

                else:
                    with self.assertRaises(SystemExit):
                        ct._update_test_status(test, ct._phases[idx+1], TEST_PEND_STATUS)

                    with self.assertRaises(SystemExit):
                        ct._update_test_status(test, phase, TEST_PASS_STATUS)

                    ct._update_test_status(test, phase, TEST_PEND_STATUS)
                    self.assertFalse(ct._is_broken(test))
                    self.assertTrue(ct._work_remains(test))

                    with self.assertRaises(SystemExit):
                        ct._update_test_status(test, phase, TEST_PEND_STATUS)

                    ct._update_test_status(test, phase, TEST_PASS_STATUS)

                    with self.assertRaises(SystemExit):
                        ct._update_test_status(test, phase, TEST_FAIL_STATUS)

                    self.assertFalse(ct._is_broken(test))
                    self.assertTrue(ct._work_remains(test))

    ###########################################################################
    def test_b_full(self):
    ###########################################################################
        tests = get_tests.get_full_test_names(["cime_test_only"], self._machine, self._compiler)
        test_id="%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        ct = TestScheduler(tests, test_id=test_id, no_batch=NO_BATCH, test_root=TEST_ROOT,
                           output_root=TEST_ROOT,compiler=self._compiler, mpilib=TEST_MPILIB)

        build_fail_test     = [item for item in tests if "TESTBUILDFAIL_" in item][0]
        build_fail_exc_test = [item for item in tests if "TESTBUILDFAILEXC" in item][0]
        run_fail_test       = [item for item in tests if "TESTRUNFAIL_" in item][0]
        run_fail_exc_test   = [item for item in tests if "TESTRUNFAILEXC" in item][0]
        pass_test           = [item for item in tests if "TESTRUNPASS" in item][0]
        test_diff_test      = [item for item in tests if "TESTTESTDIFF" in item][0]
        mem_fail_test       = [item for item in tests if "TESTMEMLEAKFAIL" in item][0]
        mem_pass_test       = [item for item in tests if "TESTMEMLEAKPASS" in item][0]
        st_arch_fail_test   = [item for item in tests if "TESTRUNSTARCFAIL" in item][0]

        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            ct.run_tests()
        finally:
            logging.getLogger().setLevel(log_lvl)

        self._wait_for_tests(test_id, expect_works=False)

        test_statuses = glob.glob("%s/*%s/TestStatus" % (self._testroot, test_id))
        self.assertEqual(len(tests), len(test_statuses))

        for test_status in test_statuses:
            ts = TestStatus(test_dir=os.path.dirname(test_status))
            test_name = ts.get_name()
            log_files = glob.glob("%s/%s*%s/TestStatus.log" % (self._testroot, test_name, test_id))
            self.assertEqual(len(log_files), 1, "Expected exactly one TestStatus.log file, found %d" % len(log_files))
            log_file = log_files[0]
            if (test_name == build_fail_test):


                assert_test_status(self, test_name, ts, MODEL_BUILD_PHASE, TEST_FAIL_STATUS)
                data = open(log_file, "r").read()
                self.assertTrue("Intentional fail for testing infrastructure" in data,
                                "Broken test did not report build error:\n%s" % data)
            elif (test_name == build_fail_exc_test):
                data = open(log_file, "r").read()
                assert_test_status(self, test_name, ts, SHAREDLIB_BUILD_PHASE, TEST_FAIL_STATUS)
                self.assertTrue("Exception from init" in data,
                                "Broken test did not report build error:\n%s" % data)
            elif (test_name == run_fail_test):
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_FAIL_STATUS)
            elif (test_name == run_fail_exc_test):
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_FAIL_STATUS)
                data = open(log_file, "r").read()
                self.assertTrue("Exception from run_phase" in data,
                                "Broken test did not report run error:\n%s" % data)
            elif (test_name == mem_fail_test):
                assert_test_status(self, test_name, ts, MEMLEAK_PHASE, TEST_FAIL_STATUS)
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_PASS_STATUS)
            elif (test_name == test_diff_test):
                assert_test_status(self, test_name, ts, "COMPARE_base_rest", TEST_FAIL_STATUS)
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_PASS_STATUS)
            elif test_name == st_arch_fail_test:
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_PASS_STATUS)
                assert_test_status(self, test_name, ts, STARCHIVE_PHASE, TEST_FAIL_STATUS)
            else:
                self.assertTrue(test_name in [pass_test, mem_pass_test])
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_PASS_STATUS)
                if (test_name == mem_pass_test):
                    assert_test_status(self, test_name, ts, MEMLEAK_PHASE, TEST_PASS_STATUS)

    ###########################################################################
    def test_c_use_existing(self):
    ###########################################################################
        tests = get_tests.get_full_test_names(["TESTBUILDFAIL_P1.f19_g16_rx1.A", "TESTRUNFAIL_P1.f19_g16_rx1.A", "TESTRUNPASS_P1.f19_g16_rx1.A"],
                                                      self._machine, self._compiler)
        test_id="%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        ct = TestScheduler(tests, test_id=test_id, no_batch=NO_BATCH, test_root=TEST_ROOT,
                           output_root=TEST_ROOT,compiler=self._compiler, mpilib=TEST_MPILIB)

        build_fail_test     = [item for item in tests if "TESTBUILDFAIL" in item][0]
        run_fail_test       = [item for item in tests if "TESTRUNFAIL" in item][0]
        pass_test           = [item for item in tests if "TESTRUNPASS" in item][0]

        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            ct.run_tests()
        finally:
            logging.getLogger().setLevel(log_lvl)

        test_statuses = glob.glob("%s/*%s/TestStatus" % (self._testroot, test_id))
        self.assertEqual(len(tests), len(test_statuses))

        self._wait_for_tests(test_id, expect_works=False)

        for test_status in test_statuses:
            casedir = os.path.dirname(test_status)
            ts = TestStatus(test_dir=casedir)
            test_name = ts.get_name()
            if test_name == build_fail_test:
                assert_test_status(self, test_name, ts, MODEL_BUILD_PHASE, TEST_FAIL_STATUS)
                with TestStatus(test_dir=casedir) as ts:
                    ts.set_status(MODEL_BUILD_PHASE, TEST_PEND_STATUS)
            elif test_name == run_fail_test:
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_FAIL_STATUS)
                with TestStatus(test_dir=casedir) as ts:
                    ts.set_status(SUBMIT_PHASE, TEST_PEND_STATUS)
            else:
                self.assertTrue(test_name == pass_test)
                assert_test_status(self, test_name, ts, MODEL_BUILD_PHASE, TEST_PASS_STATUS)
                assert_test_status(self, test_name, ts, SUBMIT_PHASE, TEST_PASS_STATUS)
                assert_test_status(self, test_name, ts, RUN_PHASE, TEST_PASS_STATUS)

        os.environ["TESTBUILDFAIL_PASS"] = "True"
        os.environ["TESTRUNFAIL_PASS"] = "True"
        ct2 = TestScheduler(tests, test_id=test_id, no_batch=NO_BATCH, use_existing=True,
                            test_root=TEST_ROOT,output_root=TEST_ROOT,compiler=self._compiler,
                            mpilib=TEST_MPILIB)

        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            ct2.run_tests()
        finally:
            logging.getLogger().setLevel(log_lvl)

        self._wait_for_tests(test_id)

        for test_status in test_statuses:
            ts = TestStatus(test_dir=os.path.dirname(test_status))
            test_name = ts.get_name()
            assert_test_status(self, test_name, ts, MODEL_BUILD_PHASE, TEST_PASS_STATUS)
            assert_test_status(self, test_name, ts, SUBMIT_PHASE, TEST_PASS_STATUS)
            assert_test_status(self, test_name, ts, RUN_PHASE, TEST_PASS_STATUS)

        del os.environ["TESTBUILDFAIL_PASS"]
        del os.environ["TESTRUNFAIL_PASS"]

        # test that passed tests are not re-run

        ct2 = TestScheduler(tests, test_id=test_id, no_batch=NO_BATCH, use_existing=True,
                            test_root=TEST_ROOT,output_root=TEST_ROOT,compiler=self._compiler,
                            mpilib=TEST_MPILIB)

        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            ct2.run_tests()
        finally:
            logging.getLogger().setLevel(log_lvl)

        self._wait_for_tests(test_id)

        for test_status in test_statuses:
            ts = TestStatus(test_dir=os.path.dirname(test_status))
            test_name = ts.get_name()
            assert_test_status(self, test_name, ts, MODEL_BUILD_PHASE, TEST_PASS_STATUS)
            assert_test_status(self, test_name, ts, SUBMIT_PHASE, TEST_PASS_STATUS)
            assert_test_status(self, test_name, ts, RUN_PHASE, TEST_PASS_STATUS)

    ###########################################################################
    def test_d_retry(self):
    ###########################################################################
        args = ["TESTBUILDFAIL_P1.f19_g16_rx1.A", "TESTRUNFAIL_P1.f19_g16_rx1.A", "TESTRUNPASS_P1.f19_g16_rx1.A", "--retry=1"]

        self._create_test(args)

###############################################################################
class P_TestJenkinsGenericJob(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def setUp(self):
    ###########################################################################
        if CIME.utils.get_model() != "e3sm":
            self.skipTest("Skipping Jenkins tests. E3SM feature")
        TestCreateTestCommon.setUp(self)

        # Need to run in a subdir in order to not have CTest clash. Name it
        # such that it should be cleaned up by the parent tearDown
        self._testdir = os.path.join(self._testroot, "jenkins_test_%s" % self._baseline_name)
        os.makedirs(self._testdir)

        # Change root to avoid clashing with other jenkins_generic_jobs
        self._jenkins_root = os.path.join(self._testdir, "J")

    ###########################################################################
    def tearDown(self):
    ###########################################################################
        TestCreateTestCommon.tearDown(self)

        if "TESTRUNDIFF_ALTERNATE" in os.environ:
            del os.environ["TESTRUNDIFF_ALTERNATE"]

    ###########################################################################
    def simple_test(self, expect_works, extra_args, build_name=None):
    ###########################################################################
        if NO_BATCH:
            extra_args += " --no-batch"

        # Need these flags to test dashboard if e3sm
        if CIME.utils.get_model() == "e3sm" and build_name is not None:
            extra_args += " -p ACME_test --submit-to-cdash --cdash-build-group=Nightly -c %s" % build_name

        run_cmd_assert_result(self, "%s/jenkins_generic_job -r %s %s -B %s" % (TOOLS_DIR, self._testdir, extra_args, self._baseline_area),
                              from_dir=self._testdir, expected_stat=(0 if expect_works else CIME.utils.TESTS_FAILED_ERR_CODE))

    ###########################################################################
    def threaded_test(self, expect_works, extra_args, build_name=None):
    ###########################################################################
        try:
            self.simple_test(expect_works, extra_args, build_name)
        except AssertionError as e:
            self._thread_error = str(e)

    ###########################################################################
    def assert_num_leftovers(self, suite):
    ###########################################################################
        num_tests_in_tiny = len(get_tests.get_test_suite(suite))

        jenkins_dirs = glob.glob("%s/*%s*/" % (self._jenkins_root, self._baseline_name.capitalize())) # case dirs
        # scratch_dirs = glob.glob("%s/*%s*/" % (self._testroot, test_id)) # blr/run dirs

        self.assertEqual(num_tests_in_tiny, len(jenkins_dirs),
                         msg="Wrong number of leftover directories in %s, expected %d, see %s" % \
                             (self._jenkins_root, num_tests_in_tiny, jenkins_dirs))

        # JGF: Can't test this at the moment due to root change flag given to jenkins_generic_job
        # self.assertEqual(num_tests_in_tiny + 1, len(scratch_dirs),
        #                  msg="Wrong number of leftover directories in %s, expected %d, see %s" % \
        #                      (self._testroot, num_tests_in_tiny, scratch_dirs))

    ###########################################################################
    def test_jenkins_generic_job(self):
    ###########################################################################
        # Generate fresh baselines so that this test is not impacted by
        # unresolved diffs
        self.simple_test(True, "-t cime_test_only_pass -g -b %s" % self._baseline_name)
        self.assert_num_leftovers("cime_test_only_pass")

        build_name = "jenkins_generic_job_pass_%s" % CIME.utils.get_timestamp()
        self.simple_test(True, "-t cime_test_only_pass -b %s" % self._baseline_name, build_name=build_name)
        self.assert_num_leftovers("cime_test_only_pass") # jenkins_generic_job should have automatically cleaned up leftovers from prior run
        assert_dashboard_has_build(self, build_name)

    ###########################################################################
    def test_jenkins_generic_job_kill(self):
    ###########################################################################
        build_name = "jenkins_generic_job_kill_%s" % CIME.utils.get_timestamp()
        run_thread = threading.Thread(target=self.threaded_test, args=(False, " -t cime_test_only_slow_pass -b master --baseline-compare=no", build_name))
        run_thread.daemon = True
        run_thread.start()

        time.sleep(120)

        kill_subprocesses(sig=signal.SIGTERM)

        run_thread.join(timeout=10)

        self.assertFalse(run_thread.isAlive(), msg="jenkins_generic_job should have finished")
        self.assertTrue(self._thread_error is None, msg="Thread had failure: %s" % self._thread_error)
        assert_dashboard_has_build(self, build_name)

    ###########################################################################
    def test_jenkins_generic_job_realistic_dash(self):
    ###########################################################################
        # The actual quality of the cdash results for this test can only
        # be inspected manually

        # Generate fresh baselines so that this test is not impacted by
        # unresolved diffs
        self.simple_test(False, "-t cime_test_all -g -b %s" % self._baseline_name)
        self.assert_num_leftovers("cime_test_all")

        # Should create a diff
        os.environ["TESTRUNDIFF_ALTERNATE"] = "True"

        # Should create a nml diff
        # Modify namelist
        fake_nl = """
 &fake_nml
   fake_item = 'fake'
   fake = .true.
/"""
        baseline_glob = glob.glob(os.path.join(self._baseline_area, self._baseline_name, "TESTRUNPASS*"))
        self.assertEqual(len(baseline_glob), 1, msg="Expected one match, got:\n%s" % "\n".join(baseline_glob))

        import stat
        for baseline_dir in baseline_glob:
            nl_path = os.path.join(baseline_dir, "CaseDocs", "datm_in")
            self.assertTrue(os.path.isfile(nl_path), msg="Missing file %s" % nl_path)

            os.chmod(nl_path, stat.S_IRUSR | stat.S_IWUSR)
            with open(nl_path, "a") as nl_file:
                nl_file.write(fake_nl)

        build_name = "jenkins_generic_job_mixed_%s" % CIME.utils.get_timestamp()
        self.simple_test(False, "-t cime_test_all -b %s" % self._baseline_name, build_name=build_name)
        self.assert_num_leftovers("cime_test_all") # jenkins_generic_job should have automatically cleaned up leftovers from prior run
        assert_dashboard_has_build(self, build_name)

###############################################################################
class M_TestCimePerformance(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_cime_case_ctrl_performance(self):
    ###########################################################################

        ts = time.time()

        num_repeat = 5
        for _ in range(num_repeat):
            self._create_test(["cime_tiny --no-build"])

        elapsed = time.time() - ts

        print("Perf test result: {:0.2f}".format(elapsed))

###############################################################################
class T_TestRunRestart(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_run_restart(self):
    ###########################################################################
        driver = CIME.utils.get_cime_default_driver()
        if driver == "mct":
            walltime="00:15:00"
        else:
            walltime="00:30:00"

        self._create_test(["--walltime "+walltime,"NODEFAIL_P1.f09_g16.X"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "{}.{}".format(CIME.utils.get_full_test_name("NODEFAIL_P1.f09_g16.X", machine=self._machine, compiler=self._compiler), self._baseline_name))
        rundir = run_cmd_no_fail("./xmlquery RUNDIR --value", from_dir=casedir)
        fail_sentinel = os.path.join(rundir, "FAIL_SENTINEL")
        self.assertTrue(os.path.exists(fail_sentinel), msg="Missing %s" % fail_sentinel)

        self.assertEqual(open(fail_sentinel, "r").read().count("FAIL"), 3)

    ###########################################################################
    def test_run_restart_too_many_fails(self):
    ###########################################################################
        driver = CIME.utils.get_cime_default_driver()
        if driver == "mct":
            walltime="00:15:00"
        else:
            walltime="00:30:00"

        self._create_test(["--walltime "+walltime,"NODEFAIL_P1.f09_g16.X"], test_id=self._baseline_name, env_changes="NODEFAIL_NUM_FAILS=5", run_errors=True)

        casedir = os.path.join(self._testroot,
                               "{}.{}".format(CIME.utils.get_full_test_name("NODEFAIL_P1.f09_g16.X", machine=self._machine, compiler=self._compiler), self._baseline_name))
        rundir = run_cmd_no_fail("./xmlquery RUNDIR --value", from_dir=casedir)
        fail_sentinel = os.path.join(rundir, "FAIL_SENTINEL")
        self.assertTrue(os.path.exists(fail_sentinel), msg="Missing %s" % fail_sentinel)

        self.assertEqual(open(fail_sentinel, "r").read().count("FAIL"), 4)

###############################################################################
class Q_TestBlessTestResults(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def tearDown(self):
    ###########################################################################
        TestCreateTestCommon.tearDown(self)

        if "TESTRUNDIFF_ALTERNATE" in os.environ:
            del os.environ["TESTRUNDIFF_ALTERNATE"]

    ###############################################################################
    def test_bless_test_results(self):
    ###############################################################################
        # Generate some baselines
        test_name = "TESTRUNDIFF_P1.f19_g16_rx1.A"

        if CIME.utils.get_model() == "e3sm":
            genargs = ["-g", "-o", "-b", self._baseline_name, test_name]
            compargs = ["-c", "-b", self._baseline_name, test_name]
        else:
            genargs = ["-g", self._baseline_name, "-o", test_name,
                       "--baseline-root ", self._baseline_area]
            compargs = ["-c", self._baseline_name, test_name,
                       "--baseline-root ", self._baseline_area]
        self._create_test(genargs)

        # Hist compare should pass
        self._create_test(compargs)

        # Change behavior
        os.environ["TESTRUNDIFF_ALTERNATE"] = "True"

        # Hist compare should now fail
        test_id = "%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        self._create_test(compargs, test_id=test_id, run_errors=True)

        # compare_test_results should detect the fail
        cpr_cmd = "{}/compare_test_results --test-root {} -t {} 2>&1" \
                  .format(TOOLS_DIR, TEST_ROOT, test_id)
        output = run_cmd_assert_result(self, cpr_cmd, expected_stat=CIME.utils.TESTS_FAILED_ERR_CODE)

        # use regex
        expected_pattern = re.compile(r'FAIL %s[^\s]* BASELINE' % test_name)
        the_match = expected_pattern.search(output)
        self.assertNotEqual(the_match, None,
                            msg="Cmd '%s' failed to display failed test in output:\n%s" % (cpr_cmd, output))

        # Bless
        run_cmd_no_fail("{}/bless_test_results --test-root {} --hist-only --force -t {}"
                        .format(TOOLS_DIR, TEST_ROOT, test_id))

        # Hist compare should now pass again
        self._create_test(compargs)

    ###############################################################################
    def test_rebless_namelist(self):
    ###############################################################################
        # Generate some namelist baselines
        test_to_change = "TESTRUNPASS_P1.f19_g16_rx1.A"
        if CIME.utils.get_model() == "e3sm":
            genargs = ["-n", "-g", "-o", "-b", self._baseline_name, "cime_test_only_pass"]
            compargs = ["-n", "-c", "-b", self._baseline_name, "cime_test_only_pass"]
        else:
            genargs = ["-n", "-g", self._baseline_name, "-o",  "cime_test_only_pass"]
            compargs = ["-n", "-c", self._baseline_name, "cime_test_only_pass"]

        self._create_test(genargs)

        # Basic namelist compare
        test_id = "%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        self._create_test(compargs, test_id=test_id)

        # Check standalone case.cmpgen_namelists
        casedir = os.path.join(self._testroot,
                               "%s.C.%s" % (CIME.utils.get_full_test_name(test_to_change, machine=self._machine, compiler=self._compiler), test_id))
        run_cmd_assert_result(self, "./case.cmpgen_namelists", from_dir=casedir)

        # compare_test_results should pass
        cpr_cmd = "{}/compare_test_results --test-root {} -n -t {} 2>&1" \
            .format(TOOLS_DIR, TEST_ROOT, test_id)
        output = run_cmd_assert_result(self, cpr_cmd)

        # use regex
        expected_pattern = re.compile(r'PASS %s[^\s]* NLCOMP' % test_to_change)
        the_match = expected_pattern.search(output)
        self.assertNotEqual(the_match, None,
                            msg="Cmd '%s' failed to display passed test in output:\n%s" % (cpr_cmd, output))


        # Modify namelist
        fake_nl = """
 &fake_nml
   fake_item = 'fake'
   fake = .true.
/"""
        baseline_area = self._baseline_area
        baseline_glob = glob.glob(os.path.join(baseline_area, self._baseline_name, "TEST*"))
        self.assertEqual(len(baseline_glob), 3, msg="Expected three matches, got:\n%s" % "\n".join(baseline_glob))

        import stat
        for baseline_dir in baseline_glob:
            nl_path = os.path.join(baseline_dir, "CaseDocs", "datm_in")
            self.assertTrue(os.path.isfile(nl_path), msg="Missing file %s" % nl_path)

            os.chmod(nl_path, stat.S_IRUSR | stat.S_IWUSR)
            with open(nl_path, "a") as nl_file:
                nl_file.write(fake_nl)

        # Basic namelist compare should now fail
        test_id = "%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        self._create_test(compargs, test_id=test_id, pre_run_errors=True)
        casedir = os.path.join(self._testroot,
                               "%s.C.%s" % (CIME.utils.get_full_test_name(test_to_change, machine=self._machine, compiler=self._compiler), test_id))
        run_cmd_assert_result(self, "./case.cmpgen_namelists", from_dir=casedir, expected_stat=100)

        # preview namelists should work
        run_cmd_assert_result(self, "./preview_namelists", from_dir=casedir)

        # This should still fail
        run_cmd_assert_result(self, "./case.cmpgen_namelists", from_dir=casedir, expected_stat=100)

        # compare_test_results should fail
        cpr_cmd = "{}/compare_test_results --test-root {} -n -t {} 2>&1" \
            .format(TOOLS_DIR, TEST_ROOT, test_id)
        output = run_cmd_assert_result(self, cpr_cmd, expected_stat=CIME.utils.TESTS_FAILED_ERR_CODE)

        # use regex
        expected_pattern = re.compile(r'FAIL %s[^\s]* NLCOMP' % test_to_change)
        the_match = expected_pattern.search(output)
        self.assertNotEqual(the_match, None,
                            msg="Cmd '%s' failed to display passed test in output:\n%s" % (cpr_cmd, output))

        # Bless
        run_cmd_no_fail("{}/bless_test_results --test-root {} -n --force -t {}"
                        .format(TOOLS_DIR, TEST_ROOT, test_id))

        # Basic namelist compare should now pass again
        self._create_test(compargs)

class X_TestQueryConfig(unittest.TestCase):
    def test_query_compsets(self):
        run_cmd_no_fail("{}/query_config --compsets".format(SCRIPT_DIR))

    def test_query_components(self):
        run_cmd_no_fail("{}/query_config --components".format(SCRIPT_DIR))

    def test_query_grids(self):
        run_cmd_no_fail("{}/query_config --grids".format(SCRIPT_DIR))

    def test_query_machines(self):
        run_cmd_no_fail("{}/query_config --machines".format(SCRIPT_DIR))

###############################################################################
class Z_FullSystemTest(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_full_system(self):
    ###########################################################################
        # Put this inside any test that's slow
        if (FAST_ONLY):
            self.skipTest("Skipping slow test")

        self._create_test(["--walltime=0:15:00", "cime_developer"], test_id=self._baseline_name)

        run_cmd_assert_result(self, "%s/cs.status.%s" % (self._testroot, self._baseline_name),
                              from_dir=self._testroot)

        # Ensure that we can get test times
        test_statuses = glob.glob(os.path.join(self._testroot, "*%s" % self._baseline_name, "TestStatus"))
        for test_status in test_statuses:
            test_time = CIME.wait_for_tests.get_test_time(os.path.dirname(test_status))
            self.assertIs(type(test_time), int, msg="get time did not return int for %s" % test_status)
            self.assertTrue(test_time > 0, msg="test time was zero for %s" % test_status)

        # Test that re-running works
        tests = get_tests.get_test_suite("cime_developer", machine=self._machine, compiler=self._compiler)
        for test in tests:
            casedir = os.path.join(TEST_ROOT, "%s.%s" % (test, self._baseline_name))

            # Subtle issue: The run phases of these tests will be in the PASS state until
            # the submitted case.test script is run, which could take a while if the system is
            # busy. This potentially leaves a window where the wait_for_tests command below will
            # not wait for the re-submitted jobs to run because it sees the original PASS.
            # The code below forces things back to PEND to avoid this race condition. Note
            # that we must use the MEMLEAK phase, not the RUN phase, because RUN being in a non-PEND
            # state is how system tests know they are being re-run and must reset certain
            # case settings.
            if self._hasbatch:
                with TestStatus(test_dir=casedir) as ts:
                    ts.set_status(MEMLEAK_PHASE, TEST_PEND_STATUS)

            run_cmd_assert_result(self, "./case.submit --skip-preview-namelist", from_dir=casedir)

        self._wait_for_tests(self._baseline_name)

###############################################################################
class K_TestCimeCase(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_cime_case(self):
    ###########################################################################
        self._create_test(["--no-build", "TESTRUNPASS_P1.f19_g16_rx1.A"], test_id=self._baseline_name)

        self.assertEqual(type(MACHINE.get_value("MAX_TASKS_PER_NODE")), int)
        self.assertTrue(type(MACHINE.get_value("PROJECT_REQUIRED")) in [type(None) , bool])

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS_P1.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        with Case(casedir, read_only=False) as case:
            build_complete = case.get_value("BUILD_COMPLETE")
            self.assertFalse(build_complete,
                             msg="Build complete had wrong value '%s'" %
                             build_complete)

            case.set_value("BUILD_COMPLETE", True)
            build_complete = case.get_value("BUILD_COMPLETE")
            self.assertTrue(build_complete,
                            msg="Build complete had wrong value '%s'" %
                            build_complete)

            case.flush()

            build_complete = run_cmd_no_fail("./xmlquery BUILD_COMPLETE -value",
                                             from_dir=casedir)
            self.assertEqual(build_complete, "TRUE",
                            msg="Build complete had wrong value '%s'" %
                            build_complete)

            # Test some test properties
            self.assertEqual(case.get_value("TESTCASE"), "TESTRUNPASS")

    def _batch_test_fixture(self, testcase_name):
        if not MACHINE.has_batch_system() or NO_BATCH:
            self.skipTest("Skipping testing user prerequisites without batch systems")
        testdir = os.path.join(TEST_ROOT, testcase_name)
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        run_cmd_assert_result(self, ("{}/create_newcase --case {} --script-root {} " +

                                     "--compset X --res f19_g16 --handle-preexisting-dirs=r --output-root {}").format(
                                         SCRIPT_DIR, testcase_name, testdir, testdir),
                              from_dir=SCRIPT_DIR)
        return testdir

    ###########################################################################
    def test_cime_case_prereq(self):
    ###########################################################################
        testcase_name = 'prereq_test'
        testdir = self._batch_test_fixture(testcase_name)
        with Case(testdir, read_only=False) as case:
            if case.get_value("depend_string") is None:
                self.skipTest("Skipping prereq test, depend_string was not provided for this batch system")
            job_name = "case.run"
            prereq_name = 'prereq_test'
            batch_commands = case.submit_jobs(prereq=prereq_name, job=job_name, skip_pnl=True, dry_run=True)
            self.assertTrue(isinstance(batch_commands, collections.Sequence), "case.submit_jobs did not return a sequence for a dry run")
            self.assertTrue(len(batch_commands) > 0, "case.submit_jobs did not return any job submission string")
            # The first element in the internal sequence should just be the job name
            # The second one (batch_cmd_index) should be the actual batch submission command
            batch_cmd_index = 1
            # The prerequisite should be applied to all jobs, though we're only expecting one
            for batch_cmd in batch_commands:
                self.assertTrue(isinstance(batch_cmd, collections.Sequence), "case.submit_jobs did not return a sequence of sequences")
                self.assertTrue(len(batch_cmd) > batch_cmd_index, "case.submit_jobs returned internal sequences with length <= {}".format(batch_cmd_index))
                self.assertTrue(isinstance(batch_cmd[1], six.string_types), "case.submit_jobs returned internal sequences without the batch command string as the second parameter: {}".format(batch_cmd[1]))
                batch_cmd_args = batch_cmd[1]

                jobid_ident = "jobid"
                dep_str_fmt = case.get_env('batch').get_value('depend_string', subgroup=None)
                self.assertTrue(jobid_ident in dep_str_fmt, "dependency string doesn't include the jobid identifier {}".format(jobid_ident))
                dep_str = dep_str_fmt[:dep_str_fmt.index(jobid_ident)]

                prereq_substr = None
                while dep_str in batch_cmd_args:
                    dep_id_pos = batch_cmd_args.find(dep_str) + len(dep_str)
                    batch_cmd_args = batch_cmd_args[dep_id_pos:]
                    prereq_substr = batch_cmd_args[:len(prereq_name)]
                    if prereq_substr == prereq_name:
                        break

                self.assertTrue(prereq_name in prereq_substr, "Dependencies added, but not the user specified one")

    ###########################################################################
    def test_cime_case_allow_failed_prereq(self):
    ###########################################################################
        testcase_name = 'allow_failed_prereq_test'
        testdir = self._batch_test_fixture(testcase_name)
        with Case(testdir, read_only=False) as case:
            depend_allow = case.get_value("depend_allow_string")
            if depend_allow is None:
                self.skipTest("Skipping allow_failed_prereq test, depend_allow_string was not provided for this batch system")
            job_name = "case.run"
            prereq_name = "prereq_allow_fail_test"
            depend_allow = depend_allow.replace("jobid", prereq_name)
            batch_commands = case.submit_jobs(prereq=prereq_name, allow_fail=True, job=job_name, skip_pnl=True, dry_run=True)
            self.assertTrue(isinstance(batch_commands, collections.Sequence), "case.submit_jobs did not return a sequence for a dry run")
            num_submissions = 1
            if case.get_value("DOUT_S"):
                num_submissions = 2
            self.assertTrue(len(batch_commands) == num_submissions, "case.submit_jobs did not return any job submission strings")
            self.assertTrue(depend_allow in batch_commands[0][1])

    ###########################################################################
    def test_cime_case_resubmit_immediate(self):
    ###########################################################################
        testcase_name = 'resubmit_immediate_test'
        testdir = self._batch_test_fixture(testcase_name)
        with Case(testdir, read_only=False) as case:
            depend_string = case.get_value("depend_string")
            if depend_string is None:
                self.skipTest("Skipping resubmit_immediate test, depend_string was not provided for this batch system")
            depend_string = depend_string.replace("jobid", "")
            job_name = "case.run"
            num_submissions = 6
            case.set_value("RESUBMIT", num_submissions - 1)
            batch_commands = case.submit_jobs(job=job_name, skip_pnl=True, dry_run=True, resubmit_immediate=True)
            self.assertTrue(isinstance(batch_commands, collections.Sequence), "case.submit_jobs did not return a sequence for a dry run")
            if case.get_value("DOUT_S"):
                num_submissions = 12
            self.assertTrue(len(batch_commands) == num_submissions, "case.submit_jobs did not return {} submitted jobs".format(num_submissions))
            for i, cmd in enumerate(batch_commands):
                if i > 0:
                    self.assertTrue(depend_string in cmd[1])

    ###########################################################################
    def test_cime_case_st_archive_resubmit(self):
    ###########################################################################
        testcase_name = "st_archive_resubmit_test"
        testdir = self._batch_test_fixture(testcase_name)
        with Case(testdir, read_only=False) as case:
            case.case_setup(clean=False, test_mode=False, reset=True)
            orig_resubmit = 2
            case.set_value("RESUBMIT", orig_resubmit)
            case.case_st_archive(resubmit=False)
            new_resubmit = case.get_value("RESUBMIT")
            self.assertTrue(orig_resubmit == new_resubmit, "st_archive resubmitted when told not to")
            case.case_st_archive(resubmit=True)
            new_resubmit = case.get_value("RESUBMIT")
            self.assertTrue((orig_resubmit - 1) == new_resubmit, "st_archive did not resubmit when told to")

    ###########################################################################
    def test_cime_case_build_threaded_1(self):
    ###########################################################################
        self._create_test(["--no-build", "TESTRUNPASS_P1x1.f19_g16_rx1.A"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS_P1x1.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        with Case(casedir, read_only=False) as case:
            build_threaded = case.get_value("SMP_PRESENT")
            self.assertFalse(build_threaded)

            build_threaded = case.get_build_threaded()
            self.assertFalse(build_threaded)

            case.set_value("FORCE_BUILD_SMP", True)

            build_threaded = case.get_build_threaded()
            self.assertTrue(build_threaded)

    ###########################################################################
    def test_cime_case_build_threaded_2(self):
    ###########################################################################
        self._create_test(["--no-build", "TESTRUNPASS_P1x2.f19_g16_rx1.A"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS_P1x2.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        with Case(casedir, read_only=False) as case:
            build_threaded = case.get_value("SMP_PRESENT")
            self.assertTrue(build_threaded)

            build_threaded = case.get_build_threaded()
            self.assertTrue(build_threaded)

    ###########################################################################
    def test_cime_case_mpi_serial(self):
    ###########################################################################
        self._create_test(["--no-build", "TESTRUNPASS_Mmpi-serial_P10.f19_g16_rx1.A"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS_Mmpi-serial_P10.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        with Case(casedir, read_only=True) as case:

            # Serial cases should not be using pnetcdf
            self.assertEqual(case.get_value("CPL_PIO_TYPENAME"), "netcdf")

            # Serial cases should be using 1 task
            self.assertEqual(case.get_value("TOTALPES"), 1)

            self.assertEqual(case.get_value("NTASKS_CPL"), 1)

    ###########################################################################
    def test_cime_case_force_pecount(self):
    ###########################################################################
        self._create_test(["--no-build", "--force-procs=16", "--force-threads=8", "TESTRUNPASS.f19_g16_rx1.A"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS_P16x8.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        with Case(casedir, read_only=True) as case:
            self.assertEqual(case.get_value("NTASKS_CPL"), 16)

            self.assertEqual(case.get_value("NTHRDS_CPL"), 8)

    ###########################################################################
    def test_cime_case_xmlchange_append(self):
    ###########################################################################
        self._create_test(["--no-build", "TESTRUNPASS_P1x1.f19_g16_rx1.A"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS_P1x1.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        run_cmd_assert_result(self, "./xmlchange --id PIO_CONFIG_OPTS --val='-opt1'", from_dir=casedir)
        result = run_cmd_assert_result(self, "./xmlquery --value PIO_CONFIG_OPTS", from_dir=casedir)
        self.assertEqual(result, "-opt1")

        run_cmd_assert_result(self, "./xmlchange --id PIO_CONFIG_OPTS --val='-opt2' --append", from_dir=casedir)
        result = run_cmd_assert_result(self, "./xmlquery --value PIO_CONFIG_OPTS", from_dir=casedir)
        self.assertEqual(result, "-opt1 -opt2")

    ###########################################################################
    def test_cime_case_test_walltime_mgmt_1(self):
    ###########################################################################
        if CIME.utils.get_model() != "e3sm":
            self.skipTest("Skipping walltime test. Depends on E3SM batch settings")

        test_name = "ERS.f19_g16_rx1.A"
        machine, compiler = "blues", "gnu"
        self._create_test(["--no-setup", "--machine={}".format(machine), test_name], test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        result = run_cmd_assert_result(self, "./xmlquery JOB_WALLCLOCK_TIME --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "0:10:00")

        result = run_cmd_assert_result(self, "./xmlquery JOB_QUEUE --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "batch")

    ###########################################################################
    def test_cime_case_test_walltime_mgmt_2(self):
    ###########################################################################
        if CIME.utils.get_model() != "e3sm":
            self.skipTest("Skipping walltime test. Depends on E3SM batch settings")

        test_name = "ERS_P64.f19_g16_rx1.A"
        machine, compiler = "blues", "gnu"
        self._create_test(["--no-setup", "--machine={}".format(machine), test_name], test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        result = run_cmd_assert_result(self, "./xmlquery JOB_WALLCLOCK_TIME --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "03:00:00")

        result = run_cmd_assert_result(self, "./xmlquery JOB_QUEUE --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "batch")

    ###########################################################################
    def test_cime_case_test_walltime_mgmt_3(self):
    ###########################################################################
        if CIME.utils.get_model() != "e3sm":
            self.skipTest("Skipping walltime test. Depends on E3SM batch settings")

        test_name = "ERS_P64.f19_g16_rx1.A"
        machine, compiler = "blues", "gnu"
        self._create_test(["--no-setup", "--machine={}".format(machine), "--walltime=0:10:00", test_name], test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        result = run_cmd_assert_result(self, "./xmlquery JOB_WALLCLOCK_TIME --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "0:10:00")

        result = run_cmd_assert_result(self, "./xmlquery JOB_QUEUE --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "batch") # Not smart enough to select faster queue

    ###########################################################################
    def test_cime_case_test_walltime_mgmt_4(self):
    ###########################################################################
        if CIME.utils.get_model() != "e3sm":
            self.skipTest("Skipping walltime test. Depends on E3SM batch settings")

        test_name = "ERS_P1.f19_g16_rx1.A"
        machine, compiler = "blues", "gnu"
        self._create_test(["--no-setup", "--machine={}".format(machine), "--walltime=2:00:00", test_name], test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        result = run_cmd_assert_result(self, "./xmlquery JOB_WALLCLOCK_TIME --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "2:00:00")

        result = run_cmd_assert_result(self, "./xmlquery JOB_QUEUE --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "batch")

    ###########################################################################
    def test_cime_case_test_walltime_mgmt_5(self):
    ###########################################################################
        if CIME.utils.get_model() != "e3sm":
            self.skipTest("Skipping walltime test. Depends on E3SM batch settings")

        test_name = "ERS_P1.f19_g16_rx1.A"
        machine, compiler = "blues", "gnu"
        self._create_test(["--no-setup", "--machine={}".format(machine), test_name], test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        run_cmd_assert_result(self, "./xmlchange JOB_QUEUE=slartibartfast --subgroup=case.test", from_dir=casedir, expected_stat=1)

        run_cmd_assert_result(self, "./xmlchange JOB_QUEUE=slartibartfast --force --subgroup=case.test", from_dir=casedir)

        result = run_cmd_assert_result(self, "./xmlquery JOB_WALLCLOCK_TIME --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "03:00:00")

        result = run_cmd_assert_result(self, "./xmlquery JOB_QUEUE --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "slartibartfast")

    ###########################################################################
    def test_cime_case_test_walltime_mgmt_6(self):
    ###########################################################################
        if not self._hasbatch:
            self.skipTest("Skipping walltime test. Depends on batch system")

        test_name = "ERS_P1.f19_g16_rx1.A"
        self._create_test(["--no-build", test_name], test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        run_cmd_assert_result(self, "./xmlchange JOB_WALLCLOCK_TIME=421:32:11 --subgroup=case.test", from_dir=casedir)

        run_cmd_assert_result(self, "./case.setup --reset", from_dir=casedir)

        result = run_cmd_assert_result(self, "./xmlquery JOB_WALLCLOCK_TIME --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "421:32:11")

    ###########################################################################
    def test_cime_case_test_walltime_mgmt_7(self):
    ###########################################################################
        if not self._hasbatch:
            self.skipTest("Skipping walltime test. Depends on batch system")

        test_name = "ERS_P1.f19_g16_rx1.A"
        self._create_test(["--no-build", "--walltime=01:00:00", test_name], test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        run_cmd_assert_result(self, "./xmlchange JOB_WALLCLOCK_TIME=421:32:11 --subgroup=case.test", from_dir=casedir)

        run_cmd_assert_result(self, "./case.setup --reset", from_dir=casedir)

        result = run_cmd_assert_result(self, "./xmlquery JOB_WALLCLOCK_TIME --subgroup=case.test --value", from_dir=casedir)
        self.assertEqual(result, "421:32:11")

    ###########################################################################
    def test_cime_case_test_custom_project(self):
    ###########################################################################
        test_name = "ERS_P1.f19_g16_rx1.A"
        machine, compiler = "melvin", "gnu" # have to use a machine both models know and one that doesn't put PROJECT in any key paths
        self._create_test(["--no-setup", "--machine={}".format(machine), "--compiler={}".format(compiler), "--project=testproj", test_name],
                          test_id=self._baseline_name,
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        result = run_cmd_assert_result(self, "./xmlquery --value PROJECT --subgroup=case.test", from_dir=casedir)
        self.assertEqual(result, "testproj")

    ###########################################################################
    def test_create_test_longname(self):
    ###########################################################################
        self._create_test(["SMS.f19_g16.2000_SATM_XLND_SICE_SOCN_XROF_XGLC_SWAV", "--no-build"])

    ###########################################################################
    def test_env_loading(self):
    ###########################################################################
        if self._machine != "melvin":
            self.skipTest("Skipping env load test - Only works on melvin")

        self._create_test(["--no-build", "TESTRUNPASS.f19_g16_rx1.A"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        with Case(casedir, read_only=True) as case:
            env_mach = case.get_env("mach_specific")
            orig_env = dict(os.environ)

            env_mach.load_env(case)
            module_env = dict(os.environ)

            os.environ.clear()
            os.environ.update(orig_env)

            env_mach.load_env(case, force_method="generic")
            generic_env = dict(os.environ)

            os.environ.clear()
            os.environ.update(orig_env)

            problems = ""
            for mkey, mval in module_env.items():
                if mkey not in generic_env:
                    if not mkey.startswith("PS") and mkey != "OLDPWD":
                        problems += "Generic missing key: {}\n".format(mkey)
                elif mval != generic_env[mkey] and mkey not in ["_", "SHLVL", "PWD"] and not mkey.endswith("()"):
                    problems += "Value mismatch for key {}: {} != {}\n".format(mkey, repr(mval), repr(generic_env[mkey]))

            for gkey in generic_env.keys():
                if gkey not in module_env:
                    problems += "Modules missing key: {}\n".format(gkey)

            self.assertEqual(problems, "", msg=problems)

    ###########################################################################
    def test_case_submit_interface(self):
    ###########################################################################
        try:
            import imp
        except ImportError:
            print("imp not found, skipping case.submit interface test")
            return
        sys.path.append(TOOLS_DIR)
        case_submit_path = os.path.join(TOOLS_DIR, "case.submit")
        submit_interface = imp.load_source("case_submit_interface", case_submit_path)
        sys.argv = ["case.submit", "--batch-args", "'random_arguments_here.%j'",
                    "--mail-type", "fail", "--mail-user", "'random_arguments_here.%j'"]
        submit_interface._main_func(None, True)

    ###########################################################################
    def test_xml_caching(self):
    ###########################################################################
        self._create_test(["--no-build", "TESTRUNPASS.f19_g16_rx1.A"], test_id=self._baseline_name)

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        active = os.path.join(casedir, "env_run.xml")
        backup = os.path.join(casedir, "env_run.xml.bak")

        safe_copy(active, backup)

        with Case(casedir, read_only=False) as case:
            env_run = EnvRun(casedir, read_only=True)
            self.assertEqual(case.get_value("RUN_TYPE"), "startup")
            case.set_value("RUN_TYPE", "branch")
            self.assertEqual(case.get_value("RUN_TYPE"), "branch")
            self.assertEqual(env_run.get_value("RUN_TYPE"), "branch")

        with Case(casedir) as case:
            self.assertEqual(case.get_value("RUN_TYPE"), "branch")

        time.sleep(0.2)
        safe_copy(backup, active)

        with Case(casedir, read_only=False) as case:
            self.assertEqual(case.get_value("RUN_TYPE"), "startup")
            case.set_value("RUN_TYPE", "branch")

        with Case(casedir, read_only=False) as case:
            self.assertEqual(case.get_value("RUN_TYPE"), "branch")
            time.sleep(0.2)
            safe_copy(backup, active)
            case.read_xml() # Manual re-sync
            self.assertEqual(case.get_value("RUN_TYPE"), "startup")
            case.set_value("RUN_TYPE", "branch")
            self.assertEqual(case.get_value("RUN_TYPE"), "branch")

        with Case(casedir) as case:
            self.assertEqual(case.get_value("RUN_TYPE"), "branch")
            time.sleep(0.2)
            safe_copy(backup, active)
            env_run = EnvRun(casedir, read_only=True)
            self.assertEqual(env_run.get_value("RUN_TYPE"), "startup")

        with Case(casedir, read_only=False) as case:
            self.assertEqual(case.get_value("RUN_TYPE"), "startup")
            case.set_value("RUN_TYPE", "branch")

        # behind the back detection
        with self.assertRaises(SystemExit):
            with Case(casedir, read_only=False) as case:
                time.sleep(0.2)
                safe_copy(backup, active)

        with Case(casedir, read_only=False) as case:
            case.set_value("RUN_TYPE", "branch")

        with self.assertRaises(SystemExit):
            with Case(casedir) as case:
                time.sleep(0.2)
                safe_copy(backup, active)

###############################################################################
class X_TestSingleSubmit(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_single_submit(self):
    ###########################################################################
        # Skip unless on a batch system and users did not select no-batch
        if (not self._hasbatch):
            self.skipTest("Skipping single submit. Not valid without batch")
        if CIME.utils.get_model() != "e3sm":
            self.skipTest("Skipping single submit. E3SM experimental feature")
        if self._machine not in ["sandiatoss3"]:
            self.skipTest("Skipping single submit. Only works on sandiatoss3")

        # Keep small enough for now that we don't have to worry about load balancing
        self._create_test(["--single-submit", "SMS_Ln9_P8.f45_g37_rx1.A", "SMS_Ln9_P8.f19_g16_rx1.A"],
                          env_changes="unset CIME_GLOBAL_WALLTIME &&")

###############################################################################
class L_TestSaveTimings(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def simple_test(self, manual_timing=False):
    ###########################################################################
        timing_flag = "" if manual_timing else "--save-timing"
        driver = CIME.utils.get_cime_default_driver()
        if driver == "mct":
            walltime="00:15:00"
        else:
            walltime="00:30:00"
        self._create_test(["SMS_Ln9_P1.f19_g16_rx1.A", timing_flag, "--walltime="+walltime], test_id=self._baseline_name)

        statuses = glob.glob("%s/*%s/TestStatus" % (self._testroot, self._baseline_name))
        self.assertEqual(len(statuses), 1, msg="Should have had exactly one match, found %s" % statuses)
        casedir = os.path.dirname(statuses[0])

        with Case(casedir, read_only=True) as case:
            lids = get_lids(case)
            timing_dir = case.get_value("SAVE_TIMING_DIR")
            casename = case.get_value("CASE")

        self.assertEqual(len(lids), 1, msg="Expected one LID, found %s" % lids)

        if manual_timing:
            run_cmd_assert_result(self, "cd %s && %s/save_provenance postrun" % (casedir, TOOLS_DIR))

        if CIME.utils.get_model() == "e3sm":
            provenance_dirs = glob.glob(os.path.join(timing_dir, "performance_archive", getpass.getuser(), casename, lids[0] + "*"))
            self.assertEqual(len(provenance_dirs), 1, msg="provenance dirs were missing")

    ###########################################################################
    def test_save_timings(self):
    ###########################################################################
        self.simple_test()

    ###########################################################################
    def test_save_timings_manual(self):
    ###########################################################################
        self.simple_test(manual_timing=True)

# Machinery for Macros generation tests.

class MockMachines(object):

    """A mock version of the Machines object to simplify testing."""

    def __init__(self, name, os_):
        """Store the name."""
        self.name = name
        self.os = os_

    def get_machine_name(self):
        """Return the name we were given."""
        return self.name

    def get_value(self, var_name):
        """Allow the operating system to be queried."""
        assert var_name == "OS", "Build asked for a value not " \
            "implemented in the testing infrastructure."
        return self.os

    def is_valid_compiler(self, _): # pylint:disable=no-self-use
        """Assume all compilers are valid."""
        return True

    def is_valid_MPIlib(self, _):
        """Assume all MPILIB settings are valid."""
        return True

# pragma pylint: disable=unused-argument
    def get_default_MPIlib(self, attributes=None):
        return "mpich2"

    def get_default_compiler(self):
        return "intel"


def get_macros(macro_maker, build_xml, build_system):
    """Generate build system ("Macros" file) output from config_compilers XML.

    Arguments:
    macro_maker - The underlying Build object.
    build_xml - A string containing the XML to operate on.
    build_system - Either "Makefile" or "CMake", depending on desired output.

    The return value is a string containing the build system output.
    """
    # Build.write_macros expects file-like objects as input, so
    # we need to wrap the strings in StringIO objects.
    xml = six.StringIO(str(build_xml))
    output = six.StringIO()
    output_format = None
    if build_system == "Makefile":
        output_format = "make"
    elif build_system == "CMake":
        output_format = "cmake"
    else:
        output_format = build_system

    macro_maker.write_macros_file(macros_file=output,
                                  output_format=output_format, xml=xml)
    return str(output.getvalue())


def _wrap_config_compilers_xml(inner_string):
    """Utility function to create a config_compilers XML string.

    Pass this function a string containing <compiler> elements, and it will add
    the necessary header/footer to the file.
    """
    _xml_template = """<?xml version="1.0" encoding="UTF-8"?>
<config_compilers>
{}
</config_compilers>
"""

    return _xml_template.format(inner_string)


class MakefileTester(object):

    """Helper class for checking Makefile output.

    Public methods:
    __init__
    query_var
    assert_variable_equals
    assert_variable_matches
    """
# Note that the following is a Makefile and the echo line must begin with a tab
    _makefile_template = """
include Macros
query:
\techo '$({})' > query.out
"""

    def __init__(self, parent, make_string):
        """Constructor for Makefile test helper class.

        Arguments:
        parent - The TestCase object that is using this item.
        make_string - Makefile contents to test.
        """
        self.parent = parent
        self.make_string = make_string

    def query_var(self, var_name, env, var):
        """Request the value of a variable in the Makefile, as a string.

        Arguments:
        var_name - Name of the variable to query.
        env - A dict containing extra environment variables to set when calling
              make.
        var - A dict containing extra make variables to set when calling make.
              (The distinction between env and var actually matters only for
               CMake, though.)
        """
        if env is None:
            env = dict()
        if var is None:
            var = dict()

        # Write the Makefile strings to temporary files.
        temp_dir = tempfile.mkdtemp()
        macros_file_name = os.path.join(temp_dir, "Macros")
        makefile_name = os.path.join(temp_dir, "Makefile")
        output_name = os.path.join(temp_dir, "query.out")

        with open(macros_file_name, "w") as macros_file:
            macros_file.write(self.make_string)
        with open(makefile_name, "w") as makefile:
            makefile.write(self._makefile_template.format(var_name))

        environment = os.environ.copy()
        environment.update(env)
        environment.update(var)
        gmake_exe = MACHINE.get_value("GMAKE")
        if gmake_exe is None:
            gmake_exe = "gmake"
        run_cmd_assert_result(self.parent, "%s query --directory=%s 2>&1" % (gmake_exe, temp_dir), env=environment)

        with open(output_name, "r") as output:
            query_result = output.read().strip()

        # Clean up the Makefiles.
        shutil.rmtree(temp_dir)

        return query_result

    def assert_variable_equals(self, var_name, value, env=None, var=None):
        """Assert that a variable in the Makefile has a given value.

        Arguments:
        var_name - Name of variable to check.
        value - The string that the variable value should be equal to.
        env - Optional. Dict of environment variables to set when calling make.
        var - Optional. Dict of make variables to set when calling make.
        """
        self.parent.assertEqual(self.query_var(var_name, env, var), value)

    def assert_variable_matches(self, var_name, regex, env=None, var=None):
        """Assert that a variable in the Makefile matches a regex.

        Arguments:
        var_name - Name of variable to check.
        regex - The regex to match.
        env - Optional. Dict of environment variables to set when calling make.
        var - Optional. Dict of make variables to set when calling make.
        """
        self.parent.assertRegexpMatches(self.query_var(var_name, env, var), regex)


class CMakeTester(object):

    """Helper class for checking CMake output.

    Public methods:
    __init__
    query_var
    assert_variable_equals
    assert_variable_matches
    """

    _cmakelists_template = """
include(./Macros.cmake)
file(WRITE query.out "${{{}}}")
"""

    def __init__(self, parent, cmake_string):
        """Constructor for CMake test helper class.

        Arguments:
        parent - The TestCase object that is using this item.
        cmake_string - CMake contents to test.
        """
        self.parent = parent
        self.cmake_string = cmake_string

    def query_var(self, var_name, env, var):
        """Request the value of a variable in Macros.cmake, as a string.

        Arguments:
        var_name - Name of the variable to query.
        env - A dict containing extra environment variables to set when calling
              cmake.
        var - A dict containing extra CMake variables to set when calling cmake.
        """
        if env is None:
            env = dict()
        if var is None:
            var = dict()

        # Write the CMake strings to temporary files.
        temp_dir = tempfile.mkdtemp()
        macros_file_name = os.path.join(temp_dir, "Macros.cmake")
        cmakelists_name = os.path.join(temp_dir, "CMakeLists.txt")
        output_name = os.path.join(temp_dir, "query.out")

        with open(macros_file_name, "w") as macros_file:
            for key in var:
                macros_file.write("set({} {})\n".format(key, var[key]))
            macros_file.write(self.cmake_string)
        with open(cmakelists_name, "w") as cmakelists:
            cmakelists.write(self._cmakelists_template.format(var_name))

        environment = os.environ.copy()
        environment.update(env)
        os_ = MACHINE.get_value("OS")
        # cmake will not work on cray systems without this flag
        if os_ == "CNL":
            cmake_args = "-DCMAKE_SYSTEM_NAME=Catamount"
        else:
            cmake_args = ""

        run_cmd_assert_result(self.parent, "cmake %s . 2>&1" % cmake_args, from_dir=temp_dir, env=environment)

        with open(output_name, "r") as output:
            query_result = output.read().strip()

        # Clean up the CMake files.
        shutil.rmtree(temp_dir)

        return query_result

    def assert_variable_equals(self, var_name, value, env=None, var=None):
        """Assert that a variable in the CMakeLists has a given value.

        Arguments:
        var_name - Name of variable to check.
        value - The string that the variable value should be equal to.
        env - Optional. Dict of environment variables to set when calling cmake.
        var - Optional. Dict of CMake variables to set when calling cmake.
        """
        self.parent.assertEqual(self.query_var(var_name, env, var), value)

    def assert_variable_matches(self, var_name, regex, env=None, var=None):
        """Assert that a variable in the CMkeLists matches a regex.

        Arguments:
        var_name - Name of variable to check.
        regex - The regex to match.
        env - Optional. Dict of environment variables to set when calling cmake.
        var - Optional. Dict of CMake variables to set when calling cmake.
        """
        self.parent.assertRegexpMatches(self.query_var(var_name, env, var), regex)


###############################################################################
class G_TestMacrosBasic(unittest.TestCase):
###############################################################################

    """Basic infrastructure tests.

    This class contains tests that do not actually depend on the output of the
    macro file conversion. This includes basic smoke testing and tests of
    error-handling in the routine.
    """

    def test_script_is_callable(self):
        """The test script can be called on valid output without dying."""
        # This is really more a smoke test of this script than anything else.
        maker = Compilers(MockMachines("mymachine", "SomeOS"), version=2.0)
        test_xml = _wrap_config_compilers_xml("<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>")
        get_macros(maker, test_xml, "Makefile")

    def test_script_rejects_bad_xml(self):
        """The macro writer rejects input that's not valid XML."""
        maker = Compilers(MockMachines("mymachine", "SomeOS"), version=2.0)
        with self.assertRaises(ParseError):
            get_macros(maker, "This is not valid XML.", "Makefile")

    def test_script_rejects_bad_build_system(self):
        """The macro writer rejects a bad build system string."""
        maker = Compilers(MockMachines("mymachine", "SomeOS"), version=2.0)
        bad_string = "argle-bargle."
        with assertRaisesRegex(self,
                SystemExit,
                "Unrecognized build system provided to write_macros: " + bad_string):
            get_macros(maker, "This string is irrelevant.", bad_string)


###############################################################################
class H_TestMakeMacros(unittest.TestCase):
###############################################################################

    """Makefile macros tests.

    This class contains tests of the Makefile output of Build.

    Aside from the usual setUp and test methods, this class has a utility method
    (xml_to_tester) that converts XML input directly to a MakefileTester object.
    """
    def setUp(self):
        self.test_os       = "SomeOS"
        self.test_machine  = "mymachine"
        self.test_compiler = MACHINE.get_default_compiler() if TEST_COMPILER is None else TEST_COMPILER
        self.test_mpilib   = MACHINE.get_default_MPIlib(attributes={"compiler":self.test_compiler}) if TEST_MPILIB is None else TEST_MPILIB

        self._maker = Compilers(MockMachines(self.test_machine, self.test_os), version=2.0)

    def xml_to_tester(self, xml_string):
        """Helper that directly converts an XML string to a MakefileTester."""
        test_xml = _wrap_config_compilers_xml(xml_string)
        return MakefileTester(self, get_macros(self._maker, test_xml, "Makefile"))

    def test_generic_item(self):
        """The macro writer can write out a single generic item."""
        xml_string = "<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>"
        tester = self.xml_to_tester(xml_string)
        tester.assert_variable_equals("SUPPORTS_CXX", "FALSE")

    def test_machine_specific_item(self):
        """The macro writer can pick out a machine-specific item."""
        xml1 = """<compiler MACH="{}"><SUPPORTS_CXX>TRUE</SUPPORTS_CXX></compiler>""".format(self.test_machine)
        xml2 = """<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")
        # Do this a second time, but with elements in the reverse order, to
        # ensure that the code is not "cheating" by taking the first match.
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")

    def test_ignore_non_match(self):
        """The macro writer ignores an entry with the wrong machine name."""
        xml1 = """<compiler MACH="bad"><SUPPORTS_CXX>TRUE</SUPPORTS_CXX></compiler>"""
        xml2 = """<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SUPPORTS_CXX", "FALSE")
        # Again, double-check that we don't just get lucky with the order.
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("SUPPORTS_CXX", "FALSE")

    def test_os_specific_item(self):
        """The macro writer can pick out an OS-specific item."""
        xml1 = """<compiler OS="{}"><SUPPORTS_CXX>TRUE</SUPPORTS_CXX></compiler>""".format(self.test_os)
        xml2 = """<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")

    def test_mach_other_compiler(self):
        """The macro writer compiler-specific logic works as expected."""
        xml1 = """<compiler COMPILER="{}"><CFLAGS><base>a b c</base></CFLAGS></compiler>""".format(self.test_compiler)
        xml2 = """<compiler MACH="{}" COMPILER="other"><CFLAGS><base>x y z</base></CFLAGS></compiler>""".format(self.test_machine)
        xml3 = """<compiler MACH="{}" COMPILER="{}"><CFLAGS><append>x y z</append></CFLAGS></compiler>""".format(self.test_machine,self.test_compiler)
        xml4 = """<compiler MACH="{}" COMPILER="{}"><CFLAGS><base>x y z</base></CFLAGS></compiler>""".format(self.test_machine,self.test_compiler)
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("CFLAGS", "a b c",env={"COMPILER":self.test_compiler})
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("CFLAGS", "a b c",env={"COMPILER":self.test_compiler})
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("CFLAGS", "a b c",env={"COMPILER":self.test_compiler})
        tester = self.xml_to_tester(xml1+xml3)
        tester.assert_variable_equals("CFLAGS", "a b c x y z",env={"COMPILER":self.test_compiler})
        tester = self.xml_to_tester(xml1+xml4)
        tester.assert_variable_equals("CFLAGS", "x y z",env={"COMPILER":self.test_compiler})
        tester = self.xml_to_tester(xml4+xml1)
        tester.assert_variable_equals("CFLAGS", "x y z",env={"COMPILER":self.test_compiler})

    def test_mach_beats_os(self):
        """The macro writer chooses machine-specific over os-specific matches."""
        xml1 = """<compiler OS="{}"><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>""".format(self.test_os)
        xml2 = """<compiler MACH="{}"><SUPPORTS_CXX>TRUE</SUPPORTS_CXX></compiler>""".format(self.test_machine)
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")

    def test_mach_and_os_beats_mach(self):
        """The macro writer chooses the most-specific match possible."""
        xml1 = """<compiler MACH="{}"><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>""".format(self.test_machine)
        xml2 = """<compiler MACH="{}" OS="{}"><SUPPORTS_CXX>TRUE</SUPPORTS_CXX></compiler>"""
        xml2 = xml2.format(self.test_machine, self.test_os)
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE")

    def test_build_time_attribute(self):
        """The macro writer writes conditionals for build-time choices."""
        xml1 = """<compiler><MPI_PATH MPILIB="mpich">/path/to/mpich</MPI_PATH></compiler>"""
        xml2 = """<compiler><MPI_PATH MPILIB="openmpi">/path/to/openmpi</MPI_PATH></compiler>"""
        xml3 = """<compiler><MPI_PATH>/path/to/default</MPI_PATH></compiler>"""
        tester = self.xml_to_tester(xml1+xml2+xml3)
        tester.assert_variable_equals("MPI_PATH", "/path/to/default")
        tester.assert_variable_equals("MPI_PATH", "/path/to/mpich", env={"MPILIB": "mpich"})
        tester.assert_variable_equals("MPI_PATH", "/path/to/openmpi", env={"MPILIB": "openmpi"})
        tester = self.xml_to_tester(xml3+xml2+xml1)
        tester.assert_variable_equals("MPI_PATH", "/path/to/default")
        tester.assert_variable_equals("MPI_PATH", "/path/to/mpich", env={"MPILIB": "mpich"})
        tester.assert_variable_equals("MPI_PATH", "/path/to/openmpi", env={"MPILIB": "openmpi"})

    def test_reject_duplicate_defaults(self):
        """The macro writer dies if given many defaults."""
        xml1 = """<compiler><MPI_PATH>/path/to/default</MPI_PATH></compiler>"""
        xml2 = """<compiler><MPI_PATH>/path/to/other_default</MPI_PATH></compiler>"""
        with assertRaisesRegex(self,
                SystemExit,
                "Variable MPI_PATH is set ambiguously in config_compilers.xml."):
            self.xml_to_tester(xml1+xml2)

    def test_reject_duplicates(self):
        """The macro writer dies if given many matches for a given configuration."""
        xml1 = """<compiler><MPI_PATH MPILIB="mpich">/path/to/mpich</MPI_PATH></compiler>"""
        xml2 = """<compiler><MPI_PATH MPILIB="mpich">/path/to/mpich2</MPI_PATH></compiler>"""
        with assertRaisesRegex(self,
                SystemExit,
                "Variable MPI_PATH is set ambiguously in config_compilers.xml."):
            self.xml_to_tester(xml1+xml2)

    def test_reject_ambiguous(self):
        """The macro writer dies if given an ambiguous set of matches."""
        xml1 = """<compiler><MPI_PATH MPILIB="mpich">/path/to/mpich</MPI_PATH></compiler>"""
        xml2 = """<compiler><MPI_PATH DEBUG="FALSE">/path/to/mpi-debug</MPI_PATH></compiler>"""
        with assertRaisesRegex(self,
                SystemExit,
                "Variable MPI_PATH is set ambiguously in config_compilers.xml."):
            self.xml_to_tester(xml1+xml2)

    def test_compiler_changeable_at_build_time(self):
        """The macro writer writes information for multiple compilers."""
        xml1 = """<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>"""
        xml2 = """<compiler COMPILER="gnu"><SUPPORTS_CXX>TRUE</SUPPORTS_CXX></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE", env={"COMPILER": "gnu"})
        tester.assert_variable_equals("SUPPORTS_CXX", "FALSE")

    def test_base_flags(self):
        """Test that we get "base" compiler flags."""
        xml1 = """<compiler><FFLAGS><base>-O2</base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-O2")

    def test_machine_specific_base_flags(self):
        """Test selection among base compiler flag sets based on machine."""
        xml1 = """<compiler><FFLAGS><base>-O2</base></FFLAGS></compiler>"""
        xml2 = """<compiler MACH="{}"><FFLAGS><base>-O3</base></FFLAGS></compiler>""".format(self.test_machine)
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("FFLAGS", "-O3")

    def test_build_time_base_flags(self):
        """Test selection of base flags based on build-time attributes."""
        xml1 = """<compiler><FFLAGS><base>-O2</base></FFLAGS></compiler>"""
        xml2 = """<compiler><FFLAGS><base DEBUG="TRUE">-O3</base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("FFLAGS", "-O2")
        tester.assert_variable_equals("FFLAGS", "-O3", env={"DEBUG": "TRUE"})

    def test_build_time_base_flags_same_parent(self):
        """Test selection of base flags in the same parent element."""
        xml1 = """<base>-O2</base>"""
        xml2 = """<base DEBUG="TRUE">-O3</base>"""
        tester = self.xml_to_tester("<compiler><FFLAGS>"+xml1+xml2+"</FFLAGS></compiler>")
        tester.assert_variable_equals("FFLAGS", "-O2")
        tester.assert_variable_equals("FFLAGS", "-O3", env={"DEBUG": "TRUE"})
        # Check for order independence here, too.
        tester = self.xml_to_tester("<compiler><FFLAGS>"+xml2+xml1+"</FFLAGS></compiler>")
        tester.assert_variable_equals("FFLAGS", "-O2")
        tester.assert_variable_equals("FFLAGS", "-O3", env={"DEBUG": "TRUE"})

    def test_append_flags(self):
        """Test appending flags to a list."""
        xml1 = """<compiler><FFLAGS><base>-delicious</base></FFLAGS></compiler>"""
        xml2 = """<compiler><FFLAGS><append>-cake</append></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("FFLAGS", "-delicious -cake")
        # Order independence, as usual.
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("FFLAGS", "-delicious -cake")

    def test_machine_specific_append_flags(self):
        """Test appending flags that are either more or less machine-specific."""
        xml1 = """<compiler><FFLAGS><append>-delicious</append></FFLAGS></compiler>"""
        xml2 = """<compiler MACH="{}"><FFLAGS><append>-cake</append></FFLAGS></compiler>""".format(self.test_machine)
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_matches("FFLAGS", "^(-delicious -cake|-cake -delicious)$")

    def test_machine_specific_base_keeps_append_flags(self):
        """Test that machine-specific base flags don't override default append flags."""
        xml1 = """<compiler><FFLAGS><append>-delicious</append></FFLAGS></compiler>"""
        xml2 = """<compiler MACH="{}"><FFLAGS><base>-cake</base></FFLAGS></compiler>""".format(self.test_machine)
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("FFLAGS", "-cake -delicious")
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("FFLAGS", "-cake -delicious")

    def test_machine_specific_base_and_append_flags(self):
        """Test that machine-specific base flags coexist with machine-specific append flags."""
        xml1 = """<compiler MACH="{}"><FFLAGS><append>-delicious</append></FFLAGS></compiler>""".format(self.test_machine)
        xml2 = """<compiler MACH="{}"><FFLAGS><base>-cake</base></FFLAGS></compiler>""".format(self.test_machine)
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("FFLAGS", "-cake -delicious")
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("FFLAGS", "-cake -delicious")

    def test_append_flags_without_base(self):
        """Test appending flags to a value set before Macros is included."""
        xml1 = """<compiler><FFLAGS><append>-cake</append></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-delicious -cake", var={"FFLAGS": "-delicious"})

    def test_build_time_append_flags(self):
        """Test build_time selection of compiler flags."""
        xml1 = """<compiler><FFLAGS><append>-cake</append></FFLAGS></compiler>"""
        xml2 = """<compiler><FFLAGS><append DEBUG="TRUE">-and-pie</append></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("FFLAGS", "-cake")
        tester.assert_variable_matches("FFLAGS", "^(-cake -and-pie|-and-pie -cake)$", env={"DEBUG": "TRUE"})

    def test_environment_variable_insertion(self):
        """Test that ENV{..} inserts environment variables."""
        # DO it again with $ENV{} style
        xml1 = """<compiler><LDFLAGS><append>-L$ENV{NETCDF} -lnetcdf</append></LDFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("LDFLAGS", "-L/path/to/netcdf -lnetcdf",
                                      env={"NETCDF": "/path/to/netcdf"})

    def test_shell_command_insertion(self):
        """Test that $SHELL insert shell command output."""
        xml1 = """<compiler><FFLAGS><base>-O$SHELL{echo 2} -fast</base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-O2 -fast")

    def test_multiple_shell_commands(self):
        """Test that more than one $SHELL element can be used."""
        xml1 = """<compiler><FFLAGS><base>-O$SHELL{echo 2} -$SHELL{echo fast}</base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-O2 -fast")

    def test_env_and_shell_command(self):
        """Test that $ENV works inside $SHELL elements."""
        xml1 = """<compiler><FFLAGS><base>-O$SHELL{echo $ENV{OPT_LEVEL}} -fast</base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-O2 -fast", env={"OPT_LEVEL": "2"})

    def test_config_variable_insertion(self):
        """Test that $VAR insert variables from config_compilers."""
        # Construct an absurd chain of references just to sure that we don't
        # pass by accident, e.g. outputting things in the right order just due
        # to good luck in a hash somewhere.
        xml1 = """<MPI_LIB_NAME>stuff-${MPI_PATH}-stuff</MPI_LIB_NAME>"""
        xml2 = """<MPI_PATH>${MPICC}</MPI_PATH>"""
        xml3 = """<MPICC>${MPICXX}</MPICC>"""
        xml4 = """<MPICXX>${MPIFC}</MPICXX>"""
        xml5 = """<MPIFC>mpicc</MPIFC>"""
        tester = self.xml_to_tester("<compiler>"+xml1+xml2+xml3+xml4+xml5+"</compiler>")
        tester.assert_variable_equals("MPI_LIB_NAME", "stuff-mpicc-stuff")

    def test_config_reject_self_references(self):
        """Test that $VAR self-references are rejected."""
        # This is a special case of the next test, which also checks circular
        # references.
        xml1 = """<MPI_LIB_NAME>${MPI_LIB_NAME}</MPI_LIB_NAME>"""
        err_msg = r".* has bad \$VAR references. Check for circular references or variables that are used in a \$VAR but not actually defined."
        with assertRaisesRegex(self,SystemExit, err_msg):
            self.xml_to_tester("<compiler>"+xml1+"</compiler>")

    def test_config_reject_cyclical_references(self):
        """Test that cyclical $VAR references are rejected."""
        xml1 = """<MPI_LIB_NAME>${MPI_PATH}</MPI_LIB_NAME>"""
        xml2 = """<MPI_PATH>${MPI_LIB_NAME}</MPI_PATH>"""
        err_msg = r".* has bad \$VAR references. Check for circular references or variables that are used in a \$VAR but not actually defined."
        with assertRaisesRegex(self,SystemExit, err_msg):
            self.xml_to_tester("<compiler>"+xml1+xml2+"</compiler>")

    def test_variable_insertion_with_machine_specific_setting(self):
        """Test that machine-specific $VAR dependencies are correct."""
        xml1 = """<compiler><MPI_LIB_NAME>something</MPI_LIB_NAME></compiler>"""
        xml2 = """<compiler MACH="{}"><MPI_LIB_NAME>$MPI_PATH</MPI_LIB_NAME></compiler>""".format(self.test_machine)
        xml3 = """<compiler><MPI_PATH>${MPI_LIB_NAME}</MPI_PATH></compiler>"""
        err_msg = r".* has bad \$VAR references. Check for circular references or variables that are used in a \$VAR but not actually defined."
        with assertRaisesRegex(self,SystemExit, err_msg):
            self.xml_to_tester(xml1+xml2+xml3)

    def test_override_with_machine_and_new_attributes(self):
        """Test that overrides with machine-specific settings with added attributes work correctly."""
        xml1 = """
<compiler COMPILER="{}">
  <SCC>icc</SCC>
  <MPICXX>mpicxx</MPICXX>
  <MPIFC>mpif90</MPIFC>
  <MPICC>mpicc</MPICC>
</compiler>""".format(self.test_compiler)
        xml2 = """
<compiler COMPILER="{}" MACH="{}">
  <MPICXX>mpifoo</MPICXX>
  <MPIFC MPILIB="{}">mpiffoo</MPIFC>
  <MPICC MPILIB="NOT_MY_MPI">mpifouc</MPICC>
</compiler>
""".format(self.test_compiler, self.test_machine, self.test_mpilib)

        tester = self.xml_to_tester(xml1+xml2)

        tester.assert_variable_equals("SCC", "icc", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})
        tester.assert_variable_equals("MPICXX", "mpifoo", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})
        tester.assert_variable_equals("MPIFC", "mpiffoo", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})
        tester.assert_variable_equals("MPICC", "mpicc", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})

        tester = self.xml_to_tester(xml2+xml1)

        tester.assert_variable_equals("SCC", "icc", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})
        tester.assert_variable_equals("MPICXX", "mpifoo", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})
        tester.assert_variable_equals("MPIFC", "mpiffoo", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})
        tester.assert_variable_equals("MPICC", "mpicc", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})

    def test_override_with_machine_and_same_attributes(self):
        """Test that machine-specific conditional overrides with the same attribute work correctly."""
        xml1 = """
<compiler COMPILER="{}">
  <MPIFC MPILIB="{}">mpifc</MPIFC>
</compiler>""".format(self.test_compiler, self.test_mpilib)
        xml2 = """
<compiler MACH="{}" COMPILER="{}">
  <MPIFC MPILIB="{}">mpif90</MPIFC>
</compiler>
""".format(self.test_machine, self.test_compiler, self.test_mpilib)

        tester = self.xml_to_tester(xml1+xml2)

        tester.assert_variable_equals("MPIFC", "mpif90", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})

        tester = self.xml_to_tester(xml2+xml1)

        tester.assert_variable_equals("MPIFC", "mpif90", env={"COMPILER":self.test_compiler, "MPILIB":self.test_mpilib})

    def test_appends_not_overriden(self):
        """Test that machine-specific base value changes don't interfere with appends."""
        xml1="""
<compiler COMPILER="{}">
 <FFLAGS>
   <base>-base1</base>
   <append DEBUG="FALSE">-debug1</append>
 </FFLAGS>
</compiler>""".format(self.test_compiler)

        xml2="""
<compiler MACH="{}" COMPILER="{}">
 <FFLAGS>
   <base>-base2</base>
   <append DEBUG="TRUE">-debug2</append>
 </FFLAGS>
</compiler>""".format(self.test_machine, self.test_compiler)

        tester = self.xml_to_tester(xml1+xml2)

        tester.assert_variable_equals("FFLAGS", "-base2", env={"COMPILER": self.test_compiler})
        tester.assert_variable_equals("FFLAGS", "-base2 -debug2", env={"COMPILER": self.test_compiler, "DEBUG": "TRUE"})
        tester.assert_variable_equals("FFLAGS", "-base2 -debug1", env={"COMPILER": self.test_compiler, "DEBUG": "FALSE"})

        tester = self.xml_to_tester(xml2+xml1)

        tester.assert_variable_equals("FFLAGS", "-base2", env={"COMPILER": self.test_compiler})
        tester.assert_variable_equals("FFLAGS", "-base2 -debug2", env={"COMPILER": self.test_compiler, "DEBUG": "TRUE"})
        tester.assert_variable_equals("FFLAGS", "-base2 -debug1", env={"COMPILER": self.test_compiler, "DEBUG": "FALSE"})

    def test_multilevel_specificity(self):
        """Check that settings with multiple levels of machine-specificity can be resolved."""
        xml1="""
<compiler>
 <MPIFC DEBUG="FALSE">mpifc</MPIFC>
</compiler>"""

        xml2="""
<compiler OS="{}">
 <MPIFC MPILIB="{}">mpif03</MPIFC>
</compiler>""".format(self.test_os, self.test_mpilib)

        xml3="""
<compiler MACH="{}">
 <MPIFC DEBUG="TRUE">mpif90</MPIFC>
</compiler>""".format(self.test_machine)

        # To verify order-independence, test every possible ordering of blocks.
        testers = []
        testers.append(self.xml_to_tester(xml1+xml2+xml3))
        testers.append(self.xml_to_tester(xml1+xml3+xml2))
        testers.append(self.xml_to_tester(xml2+xml1+xml3))
        testers.append(self.xml_to_tester(xml2+xml3+xml1))
        testers.append(self.xml_to_tester(xml3+xml1+xml2))
        testers.append(self.xml_to_tester(xml3+xml2+xml1))

        for tester in testers:
            tester.assert_variable_equals("MPIFC", "mpif90", env={"COMPILER": self.test_compiler, "MPILIB": self.test_mpilib, "DEBUG": "TRUE"})
            tester.assert_variable_equals("MPIFC", "mpif03", env={"COMPILER": self.test_compiler, "MPILIB": self.test_mpilib, "DEBUG": "FALSE"})
            tester.assert_variable_equals("MPIFC", "mpifc", env={"COMPILER": self.test_compiler, "MPILIB": "NON_MATCHING_MPI", "DEBUG": "FALSE"})

    def test_remove_dependency_issues(self):
        """Check that overridden settings don't cause inter-variable dependencies."""
        xml1="""
<compiler>
 <MPIFC>${SFC}</MPIFC>
</compiler>"""

        xml2="""
<compiler MACH="{}">""".format(self.test_machine) + """
 <SFC>${MPIFC}</SFC>
 <MPIFC>mpif90</MPIFC>
</compiler>"""

        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SFC", "mpif90")
        tester.assert_variable_equals("MPIFC", "mpif90")

        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("SFC", "mpif90")
        tester.assert_variable_equals("MPIFC", "mpif90")


###############################################################################
class I_TestCMakeMacros(H_TestMakeMacros):
###############################################################################

    """CMake macros tests.

    This class contains tests of the CMake output of Build.

    This class simply inherits all of the methods of TestMakeOutput, but changes
    the definition of xml_to_tester to create a CMakeTester instead.
    """

    def xml_to_tester(self, xml_string):
        """Helper that directly converts an XML string to a MakefileTester."""
        test_xml = _wrap_config_compilers_xml(xml_string)
        if (NO_CMAKE):
            self.skipTest("Skipping cmake test")
        else:
            return CMakeTester(self, get_macros(self._maker, test_xml, "CMake"))

###############################################################################
class S_TestManageAndQuery(unittest.TestCase):
    """Tests various scripts to manage and query xml files"""

    def _run_and_assert_query_testlist(self, extra_args=""):
        """Ensure that query_testlist runs successfully with the given extra arguments"""
        files = Files()
        testlist_drv = files.get_value("TESTS_SPEC_FILE", {"component":"drv"})

        run_cmd_assert_result(self, "{}/query_testlists --xml-testlist {} {}".format(
            SCRIPT_DIR, testlist_drv, extra_args))

    def test_query_testlists_runs(self):
        """Make sure that query_testlists runs successfully

        This simply makes sure that query_testlists doesn't generate any errors
        when it runs. This helps ensure that changes in other utilities don't
        break query_testlists.
        """
        self._run_and_assert_query_testlist(extra_args="--show-options")

    def test_query_testlists_define_testtypes_runs(self):
        """Make sure that query_testlists runs successfully with the --define-testtypes argument"""
        self._run_and_assert_query_testlist(extra_args="--define-testtypes")

    def test_query_testlists_count_runs(self):
        """Make sure that query_testlists runs successfully with the --count argument"""
        self._run_and_assert_query_testlist(extra_args="--count")

    def test_query_testlists_list_runs(self):
        """Make sure that query_testlists runs successfully with the --list argument"""
        self._run_and_assert_query_testlist(extra_args="--list categories")

###############################################################################
class B_CheckCode(unittest.TestCase):
###############################################################################
    # Tests are generated in the main loop below
    longMessage = True

    all_results = None

def make_pylint_test(pyfile, all_files):
    def test(self):
        if B_CheckCode.all_results is None:
            B_CheckCode.all_results = check_code(all_files)
        #pylint: disable=unsubscriptable-object
        result = B_CheckCode.all_results[pyfile]
        self.assertTrue(result == "", msg=result)

    return test

def check_for_pylint():
    #pylint: disable=import-error
    from distutils.spawn import find_executable
    pylint = find_executable("pylint")
    if pylint is not None:
        output = run_cmd_no_fail("pylint --version")
        pylintver = re.search(r"pylint\s+(\d+)[.](\d+)[.](\d+)", output)
        major = int(pylintver.group(1))
        minor = int(pylintver.group(2))
    if pylint is None or major < 1 or (major == 1 and minor < 5):
        print("pylint version 1.5 or newer not found, pylint tests skipped")
        return False
    return True

def write_provenance_info():
    curr_commit = get_current_commit(repo=LIB_DIR)
    logging.info("\nTesting commit %s" % curr_commit)
    cime_model = CIME.utils.get_model()
    logging.info("Using cime_model = %s" % cime_model)
    logging.info("Testing machine = %s" % MACHINE.get_machine_name())
    if TEST_COMPILER is not None:
        logging.info("Testing compiler = %s"% TEST_COMPILER)
    if TEST_MPILIB is not None:
        logging.info("Testing mpilib = %s"% TEST_MPILIB)
    logging.info("Test root: %s\n" % TEST_ROOT)

def _main_func(description):
    global MACHINE
    global NO_CMAKE
    global FAST_ONLY
    global NO_BATCH
    global TEST_COMPILER
    global TEST_MPILIB
    global TEST_ROOT
    global GLOBAL_TIMEOUT
    global NO_TEARDOWN
    config = CIME.utils.get_cime_config()

    help_str = \
"""
{0} [TEST] [TEST]
OR
{0} --help

\033[1mEXAMPLES:\033[0m
    \033[1;32m# Run the full suite \033[0m
    > {0}

    \033[1;32m# Run all code checker tests \033[0m
    > {0} B_CheckCode

    \033[1;32m# Run test test_wait_for_test_all_pass from class M_TestWaitForTests \033[0m
    > {0} M_TestWaitForTests.test_wait_for_test_all_pass
""".format(os.path.basename(sys.argv[0]))

    parser = argparse.ArgumentParser(usage=help_str,
                                     description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--fast", action="store_true",
                        help="Skip full system tests, which saves a lot of time")

    parser.add_argument("--no-batch", action="store_true",
                        help="Do not submit jobs to batch system, run locally."
                        " If false, will default to machine setting.")

    parser.add_argument("--no-cmake", action="store_true",
                        help="Do not run cmake tests")

    parser.add_argument("--no-teardown", action="store_true",
                        help="Do not delete directories left behind by testing")

    parser.add_argument("--machine",
                        help="Select a specific machine setting for cime")

    parser.add_argument("--compiler",
                        help="Select a specific compiler setting for cime")

    parser.add_argument( "--mpilib",
                        help="Select a specific compiler setting for cime")

    parser.add_argument( "--test-root",
                        help="Select a specific test root for all cases created by the testing")

    parser.add_argument("--timeout", type=int,
                        help="Select a specific timeout for all tests")

    ns, args = parser.parse_known_args()

    # Now set the sys.argv to the unittest_args (leaving sys.argv[0] alone)
    sys.argv[1:] = args

    FAST_ONLY      = ns.fast
    NO_BATCH       = ns.no_batch
    NO_CMAKE       = ns.no_cmake
    GLOBAL_TIMEOUT = ns.timeout
    NO_TEARDOWN    = ns.no_teardown

    if ns.machine is not None:
        MACHINE = Machines(machine=ns.machine)
        os.environ["CIME_MACHINE"] = ns.machine
    elif "CIME_MACHINE" in os.environ:
        mach_name = os.environ["CIME_MACHINE"]
        MACHINE = Machines(machine=mach_name)
    elif config.has_option("create_test", "MACHINE"):
        MACHINE = Machines(machine=config.get("create_test", "MACHINE"))
    elif config.has_option("main", "MACHINE"):
        MACHINE = Machines(machine=config.get("main", "MACHINE"))
    else:
        MACHINE = Machines()

    if ns.compiler is not None:
        TEST_COMPILER = ns.compiler
    elif config.has_option("create_test", "COMPILER"):
        TEST_COMPILER = config.get("create_test", "COMPILER")
    elif config.has_option("main", "COMPILER"):
        TEST_COMPILER = config.get("main", "COMPILER")

    if ns.mpilib is not None:
        TEST_MPILIB = ns.mpilib
    elif config.has_option("create_test", "MPILIB"):
        TEST_MPILIB = config.get("create_test", "MPILIB")
    elif config.has_option("main", "MPILIB"):
        TEST_MPILIB = config.get("main", "MPILIB")

    if ns.test_root is not None:
        TEST_ROOT = ns.test_root
    elif config.has_option("create_test", "TEST_ROOT"):
        TEST_ROOT = config.get("create_test", "TEST_ROOT")
    else:
        TEST_ROOT = os.path.join(MACHINE.get_value("CIME_OUTPUT_ROOT"),
                                 "scripts_regression_test.%s"% CIME.utils.get_timestamp())

    args = lambda: None # just something to set attrs on
    for log_param in ["debug", "silent", "verbose"]:
        flag = "--%s" % log_param
        if flag in sys.argv:
            sys.argv.remove(flag)
            setattr(args, log_param, True)
        else:
            setattr(args, log_param, False)

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args, None)

    write_provenance_info()

    # Find all python files in repo and create a pylint test for each
    if check_for_pylint():
        files_to_test = get_all_checkable_files()

        for file_to_test in files_to_test:
            pylint_test = make_pylint_test(file_to_test, files_to_test)
            testname = "test_pylint_%s" % file_to_test.replace("/", "_").replace(".", "_")
            expect(not hasattr(B_CheckCode, testname), "Repeat %s" % testname)
            setattr(B_CheckCode, testname, pylint_test)

    try:
        unittest.main(verbosity=2, catchbreak=True)
    except SystemExit as e:
        if e.__str__() != "False":
            print("Detected failures, leaving directory:", TEST_ROOT)
        else:
            print("All pass, removing directory:", TEST_ROOT)
            if os.path.exists(TEST_ROOT) and not NO_TEARDOWN:
                shutil.rmtree(TEST_ROOT)

        raise

if (__name__ == "__main__"):
    _main_func(__doc__)
