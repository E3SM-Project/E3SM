#!/usr/bin/env python

import io, glob, os, re, shutil, signal, sys, tempfile, \
    threading, time, logging, unittest, getpass, string

from xml.etree.ElementTree import ParseError

LIB_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(LIB_DIR)
# Remove all pyc files to ensure we're testing the right things
import subprocess
subprocess.call('/bin/rm $(find . -name "*.pyc")', shell=True, cwd=LIB_DIR)

from CIME.utils import run_cmd, run_cmd_no_fail, get_lids, get_current_commit
import update_acme_tests
import CIME.test_scheduler, CIME.wait_for_tests
from  CIME.test_scheduler import TestScheduler
from  CIME.XML.compilers import Compilers
from  CIME.XML.machines import Machines
from  CIME.XML.files import Files
from  CIME.case import Case
from  CIME.code_checker import check_code, get_all_checkable_files
from  CIME.test_status import *

SCRIPT_DIR  = CIME.utils.get_scripts_root()
TOOLS_DIR   = os.path.join(SCRIPT_DIR,"Tools")
MACHINE     = Machines()
FAST_ONLY   = False
NO_BATCH    = False
NO_CMAKE    = False
TEST_ROOT = None

os.environ["CIME_GLOBAL_WALLTIME"] = "0:05:00"

# pragma pylint: disable=protected-access
###############################################################################
def run_cmd_assert_result(test_obj, cmd, from_dir=None, expected_stat=0, env=None):
###############################################################################
    from_dir = os.getcwd() if from_dir is None else from_dir
    stat, output, errput = run_cmd(cmd, from_dir=from_dir, env=env)
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
        # from cime/utils/python
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
                ts.set_status(core_phase, status)
                break
            else:
                ts.set_status(core_phase, TEST_PASS_STATUS)

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
    # Do not test ACME dashboard if model is CESM
    if CIME.utils.get_model() == "acme":
        time.sleep(10) # Give chance for cdash to update

        wget_file = tempfile.mktemp()

        run_cmd_no_fail("wget http://my.cdash.org/index.php?project=ACME_test -O %s" % wget_file)

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

        cls._testdirs.append(testdir)

        run_cmd_assert_result(self, "%s/create_newcase --case %s --compset X --res f19_g16 --output-root %s" % (SCRIPT_DIR, testdir, cls._testroot), from_dir=SCRIPT_DIR)
        run_cmd_assert_result(self, "./case.setup", from_dir=testdir)
        run_cmd_assert_result(self, "./case.build", from_dir=testdir)

        cls._do_teardown.append(testdir)

    def test_b_user_mods(self):
        cls = self.__class__

        testdir = os.path.join(cls._testroot, 'testusermods')
        if os.path.exists(testdir):
            shutil.rmtree(testdir)

        cls._testdirs.append(testdir)

        user_mods_dir = os.path.join(CIME.utils.get_python_libs_root(), "tests", "user_mods_test1")
        run_cmd_assert_result(self, "%s/create_newcase --case %s --compset X --res f19_g16 --user-mods-dir %s --output-root %s"
                              % (SCRIPT_DIR, testdir, user_mods_dir, cls._testroot),from_dir=SCRIPT_DIR)

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
        cls._testdirs.append(testdir)

        run_cmd_assert_result(self, "%s/create_clone --clone %s --case %s --keepexe" %
                              (SCRIPT_DIR, prevtestdir, testdir),from_dir=SCRIPT_DIR)

        cls._do_teardown.append(testdir)

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
        fakeoutputroot = string.replace(cls._testroot, os.environ.get("USER"), "this_is_not_a_user")
        run_cmd_assert_result(self, "./xmlchange CIME_OUTPUT_ROOT=%s"%fakeoutputroot,
                              from_dir=prevtestdir)

        run_cmd_assert_result(self, "%s/create_clone --clone %s --case %s" %
                              (SCRIPT_DIR, prevtestdir, testdir),from_dir=SCRIPT_DIR)

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

        cls._do_teardown.append(cls._testroot)


    @classmethod
    def tearDownClass(cls):
        do_teardown = len(cls._do_teardown) > 0 and sys.exc_info() == (None, None, None)

        for tfile in cls._testdirs:
            if tfile not in cls._do_teardown:
                print "Detected failed test or user request no teardown"
                print "Leaving case directory : %s"%tfile
            elif do_teardown:
                shutil.rmtree(tfile)

###############################################################################
class M_TestWaitForTests(unittest.TestCase):
###############################################################################

    ###########################################################################
    def setUp(self):
    ###########################################################################
        self._testroot = os.path.join(TEST_ROOT,"TestWaitForTests")
        self._testdir_all_pass    = os.path.join(self._testroot, 'scripts_regression_tests.testdir_all_pass')
        self._testdir_with_fail   = os.path.join(self._testroot, 'scripts_regression_tests.testdir_with_fail')
        self._testdir_unfinished  = os.path.join(self._testroot, 'scripts_regression_tests.testdir_unfinished')
        self._testdir_unfinished2 = os.path.join(self._testroot, 'scripts_regression_tests.testdir_unfinished2')
        testdirs = [self._testdir_all_pass, self._testdir_with_fail, self._testdir_unfinished, self._testdir_unfinished2]
        for testdir in testdirs:
            if os.path.exists(testdir):
                shutil.rmtree(testdir)
            os.makedirs(testdir)

        for r in range(10):
            for testdir in testdirs:
                os.makedirs(os.path.join(testdir, str(r)))
                make_fake_teststatus(os.path.join(testdir, str(r)), "Test_%d" % r, TEST_PASS_STATUS, RUN_PHASE)

        make_fake_teststatus(os.path.join(self._testdir_with_fail,   "5"), "Test_5", TEST_FAIL_STATUS, RUN_PHASE)
        make_fake_teststatus(os.path.join(self._testdir_unfinished,  "5"), "Test_5", TEST_PEND_STATUS, RUN_PHASE)
        make_fake_teststatus(os.path.join(self._testdir_unfinished2, "5"), "Test_5", TEST_PASS_STATUS, MODEL_BUILD_PHASE)

        # Set up proxy if possible
        self._unset_proxy = setup_proxy()

        self._thread_error = None

    ###########################################################################
    def tearDown(self):
    ###########################################################################
        shutil.rmtree(self._testdir_all_pass)
        shutil.rmtree(self._testdir_with_fail)
        shutil.rmtree(self._testdir_unfinished)

        kill_subprocesses()

        if (self._unset_proxy):
            del os.environ["http_proxy"]

    ###########################################################################
    def simple_test(self, testdir, expected_results, extra_args="", build_name=None):
    ###########################################################################
        # Need these flags to test dashboard if acme
        if CIME.utils.get_model() == "acme" and build_name is not None:
            extra_args += " -b %s" % build_name

        expected_stat = 0 if expected_results == ["PASS"]*len(expected_results) else CIME.utils.TESTS_FAILED_ERR_CODE
        output = run_cmd_assert_result(self, "%s/wait_for_tests -p ACME_test */TestStatus %s" % (TOOLS_DIR, extra_args),
                                       from_dir=testdir, expected_stat=expected_stat)

        lines = [line for line in output.splitlines() if line.startswith("Test '")]
        self.assertEqual(len(lines), 10)
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

        if CIME.utils.get_model() == "acme":
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
            self.assertTrue(r'<Test Status="failed"><Name>Test_5</Name>' in xml_contents)

            # TODO: Any further checking of xml output worth doing?

###############################################################################
class TestCreateTestCommon(unittest.TestCase):
###############################################################################

    ###########################################################################
    def setUp(self):
    ###########################################################################
        self._thread_error      = None
        self._unset_proxy       = setup_proxy()
        self._baseline_name     = "fake_testing_only_%s" % CIME.utils.get_timestamp()
        self._machine           = MACHINE.get_machine_name()
        self._baseline_area     = MACHINE.get_value("BASELINE_ROOT")
        self._testroot          = TEST_ROOT
        self._compiler          = MACHINE.get_default_compiler()
        self._hasbatch          = MACHINE.has_batch_system() and not NO_BATCH
        self._do_teardown       = True # Will never do teardown if test failed

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
            print "Detected failed test or user request no teardown"
            print "Leaving files:"
            for file_to_clean in files_to_clean:
                print " ", file_to_clean
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

###############################################################################
class O_TestTestScheduler(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_a_phases(self):
    ###########################################################################
        # exclude the MEMLEAK tests here.
        tests = update_acme_tests.get_full_test_names(["cime_test_only",
                                                       "^TESTMEMLEAKFAIL_P1.f19_g16.X",
                                                       "^TESTMEMLEAKPASS_P1.f19_g16.X",
                                                       "^TESTTESTDIFF_P1.f19_g16_rx1.A",
                                                       "^TESTBUILDFAILEXC_P1.f19_g16_rx1.A",
                                                       "^TESTRUNFAILEXC_P1.f19_g16_rx1.A"],
                                                      self._machine, self._compiler)
        self.assertEqual(len(tests), 3)
        ct = TestScheduler(tests, test_root=TEST_ROOT, output_root=TEST_ROOT)

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
                elif (phase == CIME.test_scheduler.MODEL_BUILD_PHASE):
                    ct._update_test_status(test, phase, TEST_PEND_STATUS)

                    if (test == build_fail_test):
                        ct._update_test_status(test, phase, TEST_FAIL_STATUS)
                        self.assertTrue(ct._is_broken(test))
                        self.assertFalse(ct._work_remains(test))
                    else:
                        ct._update_test_status(test, phase, TEST_PASS_STATUS)
                        self.assertFalse(ct._is_broken(test))
                        self.assertTrue(ct._work_remains(test))

                elif (phase == CIME.test_scheduler.RUN_PHASE):
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
        tests = update_acme_tests.get_full_test_names(["cime_test_only"], self._machine, self._compiler)
        test_id="%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        ct = TestScheduler(tests, test_id=test_id, no_batch=NO_BATCH, test_root=TEST_ROOT,output_root=TEST_ROOT)

        build_fail_test     = [item for item in tests if "TESTBUILDFAIL_" in item][0]
        build_fail_exc_test = [item for item in tests if "TESTBUILDFAILEXC" in item][0]
        run_fail_test       = [item for item in tests if "TESTRUNFAIL_" in item][0]
        run_fail_exc_test   = [item for item in tests if "TESTRUNFAILEXC" in item][0]
        pass_test           = [item for item in tests if "TESTRUNPASS" in item][0]
        test_diff_test      = [item for item in tests if "TESTTESTDIFF" in item][0]
        mem_fail_test       = [item for item in tests if "TESTMEMLEAKFAIL" in item][0]
        mem_pass_test       = [item for item in tests if "TESTMEMLEAKPASS" in item][0]

        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            ct.run_tests()
        finally:
            logging.getLogger().setLevel(log_lvl)

        if (self._hasbatch):
            run_cmd_assert_result(self, "%s/wait_for_tests *%s/TestStatus" % (TOOLS_DIR, test_id), from_dir=self._testroot,
                                  expected_stat=CIME.utils.TESTS_FAILED_ERR_CODE)

        test_statuses = glob.glob("%s/*%s/TestStatus" % (self._testroot, test_id))
        self.assertEqual(len(tests), len(test_statuses))

        for test_status in test_statuses:
            ts = TestStatus(test_dir=os.path.dirname(test_status))
            test_name = ts.get_name()
            log_files = glob.glob("%s/%s*%s/TestStatus.log" % (self._testroot, test_name, test_id))
            self.assertEqual(len(log_files), 1, "Expected exactly one TestStatus.log file, found %d" % len(log_files))
            log_file = log_files[0]
            if (test_name == build_fail_test):
                self.assertEqual(ts.get_status(CIME.test_scheduler.MODEL_BUILD_PHASE), TEST_FAIL_STATUS)
                data = open(log_file, "r").read()
                self.assertTrue("Intentional fail for testing infrastructure" in data,
                                "Broken test did not report build error:\n%s" % data)
            elif (test_name == build_fail_exc_test):
                data = open(log_file, "r").read()
                self.assertEqual(ts.get_status(CIME.test_scheduler.SHAREDLIB_BUILD_PHASE), TEST_FAIL_STATUS)
                self.assertTrue("Exception from init" in data,
                                "Broken test did not report build error:\n%s" % data)
            elif (test_name == run_fail_test):
                self.assertEqual(ts.get_status(CIME.test_scheduler.RUN_PHASE), TEST_FAIL_STATUS)
            elif (test_name == run_fail_exc_test):
                self.assertEqual(ts.get_status(CIME.test_scheduler.RUN_PHASE), TEST_FAIL_STATUS)
                data = open(log_file, "r").read()
                self.assertTrue("Exception from run_phase" in data,
                                "Broken test did not report run error:\n%s" % data)
            elif (test_name == mem_fail_test):
                self.assertEqual(ts.get_status(MEMLEAK_PHASE), TEST_FAIL_STATUS)
                self.assertEqual(ts.get_status(CIME.test_scheduler.RUN_PHASE), TEST_PASS_STATUS)
            elif (test_name == test_diff_test):
                self.assertEqual(ts.get_status("COMPARE_base_rest"), TEST_FAIL_STATUS)
                self.assertEqual(ts.get_status(CIME.test_scheduler.RUN_PHASE), TEST_PASS_STATUS)
            else:
                self.assertTrue(test_name in [pass_test, mem_pass_test])
                self.assertEqual(ts.get_status(CIME.test_scheduler.RUN_PHASE), TEST_PASS_STATUS)
                if (test_name == mem_pass_test):
                    self.assertEqual(ts.get_status(MEMLEAK_PHASE), TEST_PASS_STATUS)

    ###########################################################################
    def test_c_use_existing(self):
    ###########################################################################
        tests = update_acme_tests.get_full_test_names(["TESTBUILDFAIL_P1.f19_g16_rx1.A", "TESTRUNPASS_P1.f19_g16_rx1.A"],
                                                      self._machine, self._compiler)
        test_id="%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        ct = TestScheduler(tests, test_id=test_id, no_batch=NO_BATCH, no_run=True,test_root=TEST_ROOT,output_root=TEST_ROOT)

        build_fail_test     = [item for item in tests if "TESTBUILDFAIL" in item][0]
        pass_test           = [item for item in tests if "TESTRUNPASS" in item][0]

        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            ct.run_tests()
        finally:
            logging.getLogger().setLevel(log_lvl)

        test_statuses = glob.glob("%s/*%s/TestStatus" % (self._testroot, test_id))
        self.assertEqual(len(tests), len(test_statuses))

        for test_status in test_statuses:
            ts = TestStatus(test_dir=os.path.dirname(test_status))
            test_name = ts.get_name()
            if test_name == build_fail_test:
                self.assertEqual(ts.get_status(CIME.test_scheduler.MODEL_BUILD_PHASE), TEST_FAIL_STATUS)
            else:
                self.assertTrue(test_name == pass_test)
                self.assertEqual(ts.get_status(CIME.test_scheduler.MODEL_BUILD_PHASE), TEST_PASS_STATUS)
                self.assertEqual(ts.get_status(CIME.test_scheduler.RUN_PHASE), TEST_PEND_STATUS)

        os.environ["TESTBUILDFAIL_PASS"] = "True"
        ct2 = TestScheduler(tests, test_id=test_id, no_batch=NO_BATCH, use_existing=True,
                            test_root=TEST_ROOT,output_root=TEST_ROOT)

        log_lvl = logging.getLogger().getEffectiveLevel()
        logging.disable(logging.CRITICAL)
        try:
            ct2.run_tests()
        finally:
            logging.getLogger().setLevel(log_lvl)

        if (self._hasbatch):
            run_cmd_assert_result(self, "%s/wait_for_tests *%s/TestStatus" % (TOOLS_DIR, test_id), from_dir=self._testroot)

        for test_status in test_statuses:
            ts = TestStatus(test_dir=os.path.dirname(test_status))
            self.assertEqual(ts.get_status(CIME.test_scheduler.MODEL_BUILD_PHASE), TEST_PASS_STATUS)
            self.assertEqual(ts.get_status(CIME.test_scheduler.RUN_PHASE),         TEST_PASS_STATUS)

###############################################################################
class P_TestJenkinsGenericJob(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def setUp(self):
    ###########################################################################
        if CIME.utils.get_model() != "acme":
            self.skipTest("Skipping Jenkins tests. ACME feature")
        TestCreateTestCommon.setUp(self)

        # Need to run in a subdir in order to not have CTest clash. Name it
        # such that it should be cleaned up by the parent tearDown
        self._testdir = os.path.join(self._testroot, "jenkins_generic_test_subdir_%s" % self._baseline_name)
        os.makedirs(self._testdir)

        # Change root to avoid clashing with other jenkins_generic_jobs
        self._jenkins_root = os.path.join(self._testdir, "jenkins")

    ###########################################################################
    def simple_test(self, expect_works, extra_args, build_name=None):
    ###########################################################################
        if NO_BATCH:
            extra_args += " --no-batch"

        # Need these flags to test dashboard if acme
        if CIME.utils.get_model() == "acme" and build_name is not None:
            extra_args += " -p ACME_test --submit-to-cdash --cdash-build-group=Nightly -c %s" % build_name

        run_cmd_assert_result(self, "%s/jenkins_generic_job -r %s %s" % (TOOLS_DIR, self._testdir, extra_args),
                              from_dir=self._testdir, expected_stat=(0 if expect_works else CIME.utils.TESTS_FAILED_ERR_CODE))

    ###########################################################################
    def threaded_test(self, expect_works, extra_args, build_name=None):
    ###########################################################################
        try:
            self.simple_test(expect_works, extra_args, build_name)
        except AssertionError as e:
            self._thread_error = str(e)

    ###########################################################################
    def assert_num_leftovers(self, test_id=None):
    ###########################################################################
        # There should only be two directories matching the test_id in both
        # the testroot (bld/run dump area) and jenkins root
        if (test_id is None):
            test_id = self._baseline_name
        num_tests_in_tiny = len(update_acme_tests.get_test_suite("cime_test_only_pass"))

        jenkins_dirs = glob.glob("%s/*%s*/" % (self._jenkins_root, test_id)) # case dirs
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
        # Unfortunately, this test is very long-running

        # Generate fresh baselines so that this test is not impacted by
        # unresolved diffs
        self.simple_test(True, "-t cime_test_only_pass -g -b %s" % self._baseline_name)
        self.assert_num_leftovers()

        build_name = "jenkins_generic_job_pass_%s" % CIME.utils.get_timestamp()
        self.simple_test(True, "-t cime_test_only_pass -b %s" % self._baseline_name, build_name=build_name)
        self.assert_num_leftovers() # jenkins_generic_job should have automatically cleaned up leftovers from prior run
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

###############################################################################
class Q_TestBlessTestResults(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def tearDown(self):
    ###########################################################################
        TestCreateTestCommon.tearDown(self)

        if "TESTRUNDIFF_ALTERNATE" in os.environ:
            del os.environ["TESTRUNDIFF_ALTERNATE"]

    ###########################################################################
    def simple_test(self, expect_works, extra_args):
    ###########################################################################
        if NO_BATCH:
            extra_args += " --no-batch"

        if " -n " in extra_args:
            run_cmd_assert_result(self, "%s/create_test --test-root %s --output-root %s %s"
                                  % (SCRIPT_DIR, TEST_ROOT, TEST_ROOT, extra_args),
                                  expected_stat=(0 if expect_works else CIME.utils.TESTS_FAILED_ERR_CODE))
        else:
            run_cmd_assert_result(self, "%s/create_test --test-root %s --output-root %s %s"
                                  % (SCRIPT_DIR, TEST_ROOT, TEST_ROOT, extra_args),
                                  expected_stat=(0 if expect_works or self._hasbatch else CIME.utils.TESTS_FAILED_ERR_CODE))

            if self._hasbatch:
                test_id = extra_args.split()[extra_args.split().index("-t") + 1]
                run_cmd_assert_result(self, "%s/wait_for_tests *%s/TestStatus" % (TOOLS_DIR, test_id),
                                      from_dir=self._testroot, expected_stat=(0 if expect_works else CIME.utils.TESTS_FAILED_ERR_CODE))

    ###############################################################################
    def test_bless_test_results(self):
    ###############################################################################
        # Generate some baselines
        test_name = "TESTRUNDIFF_P1.f19_g16_rx1.A"

        if CIME.utils.get_model() == "acme":
            genarg = "-g -o -b %s %s" % (self._baseline_name, test_name)
            comparg = "-c -b %s %s" % (self._baseline_name, test_name)
        else:
            genarg = "-g %s -o %s" % (self._baseline_name, test_name)
            comparg = "-c %s %s" % (self._baseline_name, test_name)

        self.simple_test(True, "%s -t %s-%s" % (genarg, self._baseline_name, CIME.utils.get_timestamp()))

        # Hist compare should pass
        self.simple_test(True, "%s -t %s-%s" % (comparg, self._baseline_name, CIME.utils.get_timestamp()))

        # Change behavior
        os.environ["TESTRUNDIFF_ALTERNATE"] = "True"

        # Hist compare should now fail
        test_id = "%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        self.simple_test(False, "%s -t %s" % (comparg, test_id))

        # compare_test_results should detect the fail
        cpr_cmd = "%s/compare_test_results --test-root %s -b %s -t %s 2>&1" \
            % (TOOLS_DIR, TEST_ROOT, self._baseline_name, test_id)
        output = run_cmd_assert_result(self, cpr_cmd, expected_stat=CIME.utils.TESTS_FAILED_ERR_CODE)

        # use regex
        expected_pattern = re.compile(r'FAIL %s[^\s]* BASELINE' % test_name)
        the_match = expected_pattern.search(output)
        self.assertNotEqual(the_match, None,
                            msg="Cmd '%s' failed to display failed test in output:\n%s" % (cpr_cmd, output))

        # Bless
        run_cmd_no_fail("%s/bless_test_results --test-root %s --hist-only --force -b %s -t %s"\
                            % (TOOLS_DIR, TEST_ROOT, self._baseline_name, test_id))

        # Hist compare should now pass again
        self.simple_test(True, "%s -t %s-%s" % (comparg, self._baseline_name, CIME.utils.get_timestamp()))

    ###############################################################################
    def test_rebless_namelist(self):
    ###############################################################################
        # Generate some namelist baselines
        test_to_change = "TESTRUNPASS_P1.f19_g16_rx1.A"
        if CIME.utils.get_model() == "acme":
            genarg = "-g -o -b %s cime_test_only_pass" % self._baseline_name
            comparg = "-c -b %s cime_test_only_pass" % self._baseline_name
        else:
            genarg = "-g %s -o cime_test_only_pass" % self._baseline_name
            comparg = "-c %s cime_test_only_pass" % self._baseline_name

        self.simple_test(True, "%s -n -t %s-%s" % (genarg, self._baseline_name, CIME.utils.get_timestamp()))

        # Basic namelist compare
        test_id = "%s-%s" % (self._baseline_name, CIME.utils.get_timestamp())
        self.simple_test(True, "%s -n -t %s" % (comparg, test_id))

        # Check standalone case.cmpgen_namelists
        casedir = os.path.join(self._testroot,
                               "%s.C.%s" % (CIME.utils.get_full_test_name(test_to_change, machine=self._machine, compiler=self._compiler), test_id))
        run_cmd_assert_result(self, "./case.cmpgen_namelists", from_dir=casedir)

        # compare_test_results should pass
        cpr_cmd = "%s/compare_test_results --test-root %s -n -b %s -t %s 2>&1" \
            % (TOOLS_DIR, TEST_ROOT, self._baseline_name, test_id)
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
        self.simple_test(False, "%s -n -t %s" % (comparg, test_id))
        casedir = os.path.join(self._testroot,
                               "%s.C.%s" % (CIME.utils.get_full_test_name(test_to_change, machine=self._machine, compiler=self._compiler), test_id))
        run_cmd_assert_result(self, "./case.cmpgen_namelists", from_dir=casedir, expected_stat=100)

        # preview namelists should work
        run_cmd_assert_result(self, "./preview_namelists", from_dir=casedir)

        # This should still fail
        run_cmd_assert_result(self, "./case.cmpgen_namelists", from_dir=casedir, expected_stat=100)

        # compare_test_results should fail
        cpr_cmd = "%s/compare_test_results --test-root %s -n -b %s -t %s 2>&1" \
            % (TOOLS_DIR, TEST_ROOT, self._baseline_name, test_id)
        output = run_cmd_assert_result(self, cpr_cmd, expected_stat=CIME.utils.TESTS_FAILED_ERR_CODE)

        # use regex
        expected_pattern = re.compile(r'FAIL %s[^\s]* NLCOMP' % test_to_change)
        the_match = expected_pattern.search(output)
        self.assertNotEqual(the_match, None,
                            msg="Cmd '%s' failed to display passed test in output:\n%s" % (cpr_cmd, output))

        # Bless
        run_cmd_no_fail("%s/bless_test_results --test-root %s -n --force -b %s -t %s" \
                            % (TOOLS_DIR, TEST_ROOT, self._baseline_name, test_id))

        # Basic namelist compare should now pass again
        self.simple_test(True, "%s -n -t %s-%s" % (comparg, self._baseline_name, CIME.utils.get_timestamp()))

###############################################################################
@unittest.skip("Disabling this test until we figure out how to integrate ACME tests and CIME xml files.")
class R_TestUpdateACMETests(unittest.TestCase):
###############################################################################

    ###########################################################################
    def setUp(self):
    ###########################################################################
        # Grab all active tests
        self._testlist_allactive = os.path.join(CIME.utils.get_model_config_root(), "allactive", "testlist_allactive.xml")
        shutil.copy2(self._testlist_allactive, ".")

    ###########################################################################
    def tearDown(self):
    ###########################################################################
        shutil.copy2("testlist_allactive.xml", self._testlist_allactive)

    ###########################################################################
    def test_update_acme_tests(self):
    ###########################################################################
        # Add some testable stuff to acme tests
        pass
        # update_acme_tests._TEST_SUITES["cime_tiny"] = \
        #     (None, (("ERS.f19_g16_rx1.A", "jgftestmodtest/test_mod"),
        #             ("NCK.f19_g16_rx1.A", "jgftestmodtest/test_mod"))
        #      )

        # try:
        #     update_acme_tests.update_acme_tests(os.path.basename(self._testlist_allactive), update_acme_tests.get_test_suites())
        # except:
        #     traceback.print_tb(sys.exc_info()[2])
        #     self.assertTrue(False, str(sys.exc_info()[1]))

        # stat = run_cmd("grep 'jgftestmodtest/test_mod' %s" % os.path.basename(self._testlist_allactive))[0]
        # self.assertEqual(stat, 0, msg="update_acme_tests did not update XML")

    ###########################################################################
    def test_update_acme_tests_test_mods(self):
    ###########################################################################
        pass
        # machine = "melvin"
        # not_my_machine = "%s_jgftest" % machine

        # # Add some testable stuff to acme tests
        # update_acme_tests._TEST_SUITES["cime_tiny"] = \
        #     (None, (("ERS.f19_g16_rx1.A", "test_mod"),
        #             ("ERS.f19_g16_rx1.B", "test_mod", machine),
        #             ("ERS.f19_g16_rx1.C", "test_mod", (machine, not_my_machine)),
        #             ("ERS.f19_g16_rx1.D", "test_mod", not_my_machine),
        #             "ERS.f19_g16_rx1.E")
        #      )

        # tests = update_acme_tests.get_test_suite("cime_tiny", compiler="gnu")

        # self.assertEqual(5, len(tests))
        # self.assertTrue("ERS.f19_g16_rx1.A.melvin_gnu.test_mod" in tests)
        # self.assertTrue("ERS.f19_g16_rx1.B.melvin_gnu.test_mod" in tests)
        # self.assertTrue("ERS.f19_g16_rx1.C.melvin_gnu.test_mod" in tests)
        # self.assertTrue("ERS.f19_g16_rx1.D.melvin_gnu" in tests)
        # self.assertTrue("ERS.f19_g16_rx1.E.melvin_gnu" in tests)

        # if (CIME.utils.does_machine_have_batch()):
        #     stat, output, errput = run_cmd("%s/wait_for_tests *%s*/TestStatus" % (TOOLS_DIR, self._baseline_name), from_dir=self._testroot)
        #     self.assertEqual(stat, 0,
        #                      msg="COMMAND SHOULD HAVE WORKED\nwait_for_tests output:\n%s\n\nerrput:\n%s\n\ncode: %d" % (output, errput, stat))


        # stat, output, errput = run_cmd("%s/cs.status.%s" % (self._testroot, self._baseline_name))
        # self.assertEqual(stat, 0,
        #                  msg="COMMAND SHOULD HAVE WORKED\ncs.status output:\n%s\n\nerrput:\n%s\n\ncode: %d" % (output, errput, stat))

###############################################################################
class Z_FullSystemTest(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_full_system(self):
    ###########################################################################
        # Put this inside any test that's slow
        if (FAST_ONLY):
            self.skipTest("Skipping slow test")

        create_test_cmd =  "%s/create_test cime_developer --test-root %s --output-root %s --walltime 0:15:00 -t %s" \
            % (SCRIPT_DIR, TEST_ROOT, TEST_ROOT, self._baseline_name)
        if NO_BATCH:
            create_test_cmd += " --no-batch"

        run_cmd_assert_result(self, create_test_cmd)

        if (self._hasbatch):
            run_cmd_assert_result(self, "%s/wait_for_tests *%s/TestStatus" % (TOOLS_DIR, self._baseline_name),
                                  from_dir=self._testroot)

        run_cmd_assert_result(self, "%s/cs.status.%s" % (self._testroot, self._baseline_name),
                              from_dir=self._testroot)

        # Ensure that we can get test times
        test_statuses = glob.glob(os.path.join(self._testroot, "*%s" % self._baseline_name, "TestStatus"))
        for test_status in test_statuses:
            test_time = CIME.wait_for_tests.get_test_time(os.path.dirname(test_status))
            self.assertIs(type(test_time), int, msg="get time did not return int for %s" % test_status)
            self.assertTrue(test_time > 0, msg="test time was zero for %s" % test_status)

###############################################################################
class K_TestCimeCase(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_cime_case(self):
    ###########################################################################
        run_cmd_assert_result(self, "%s/create_test cime_test_only -t %s --no-build --test-root %s --output-root %s"
                              % (SCRIPT_DIR, self._baseline_name, TEST_ROOT, TEST_ROOT))

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

            model_specific_val = case.get_value("CPL_SEQ_OPTION")
            if CIME.utils.get_model() == 'cesm':
                self.assertEqual(model_specific_val, "RASM_OPTION1")
            else:
                self.assertEqual(model_specific_val, "CESM1_MOD")

            # Test some test properties
            self.assertEqual(case.get_value("TESTCASE"), "TESTRUNPASS")

    ###########################################################################
    def test_cime_case_mpi_serial(self):
    ###########################################################################
        run_cmd_assert_result(self, "%s/create_test TESTRUNPASS_Mmpi-serial.f19_g16_rx1.A -t %s --no-build --test-root %s --output-root %s"
                              % (SCRIPT_DIR, self._baseline_name, self._testroot, self._testroot))

        casedir = os.path.join(self._testroot,
                               "%s.%s" % (CIME.utils.get_full_test_name("TESTRUNPASS_Mmpi-serial.f19_g16_rx1.A", machine=self._machine, compiler=self._compiler), self._baseline_name))
        self.assertTrue(os.path.isdir(casedir), msg="Missing casedir '%s'" % casedir)

        with Case(casedir, read_only=True) as case:

            # Serial cases should not be using pnetcdf
            self.assertEqual(case.get_value("CPL_PIO_TYPENAME"), "netcdf")

            # Serial cases should be using 1 task
            self.assertEqual(case.get_value("TOTALPES"), 1)

###############################################################################
class X_TestSingleSubmit(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def test_single_submit(self):
    ###########################################################################
        # Skip unless on a batch system and users did not select no-batch
        if (not self._hasbatch):
            self.skipTest("Skipping single submit. Not valid without batch")
        if CIME.utils.get_model() != "acme":
            self.skipTest("Skipping single submit. ACME experimental feature")
        if self._machine != "skybridge":
            self.skipTest("Skipping single submit. Only works on skybridge")

        # Keep small enough for now that we don't have to worry about load balancing
        run_cmd_assert_result(self,
                              "unset CIME_GLOBAL_WALLTIME && %s/create_test SMS_Ln9_P8.f45_g37_rx1.A SMS_Ln9_P8.f19_g16_rx1.A  -t %s --single-submit --test-root %s --output-root %s "
                              % (SCRIPT_DIR, self._baseline_name, TEST_ROOT, TEST_ROOT))
        run_cmd_assert_result(self, "%s/wait_for_tests *%s/TestStatus" % (TOOLS_DIR, self._baseline_name),
                              from_dir=self._testroot)

###############################################################################
class L_TestSaveTimings(TestCreateTestCommon):
###############################################################################

    ###########################################################################
    def simple_test(self, manual_timing=False):
    ###########################################################################
        timing_flag = "" if manual_timing else "--save-timing"
        create_test_cmd =  "%s/create_test SMS_Ln9_P1.f19_g16_rx1.A %s --walltime 0:15:00 -t %s --test-root %s --output-root %s " \
            % (SCRIPT_DIR, timing_flag, self._baseline_name, TEST_ROOT, TEST_ROOT)
        if NO_BATCH:
            create_test_cmd += " --no-batch"

        run_cmd_assert_result(self, create_test_cmd)
        if (self._hasbatch):
            run_cmd_assert_result(self, "%s/wait_for_tests *%s/TestStatus" % (TOOLS_DIR, self._baseline_name),
                                  from_dir=self._testroot)

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

        if CIME.utils.get_model() == "acme":
            provenance_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), casename, lids[0])
            self.assertTrue(os.path.isdir(provenance_dir), msg="'%s' was missing" % provenance_dir)

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

    def get_default_MPIlib(self):
        return "mpich2"

    def get_default_compiler(self):
        return "intel"


def get_macros(macro_maker, build_xml, build_system):
    """Generate build system ("Macros" file) output from config_build XML.

    Arguments:
    macro_maker - The underlying Build object.
    build_xml - A string containing the XML to operate on.
    build_system - Either "Makefile" or "CMake", depending on desired output.

    The return value is a string containing the build system output.
    """
    # Build.write_macros expects file-like objects as input, so
    # we need to wrap the strings in StringIO objects.
    xml = io.StringIO(unicode(build_xml))
    output = io.StringIO()
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


def _wrap_config_build_xml(inner_string):
    """Utility function to create a config_build XML string.

    Pass this function a string containing <compiler> elements, and it will add
    the necessary header/footer to the file.
    """
    _xml_template = """<?xml version="1.0" encoding="UTF-8"?>
<config_build>
{}
</config_build>
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
                macros_file.write("set(CIME_{} {})\n".format(key, var[key]))
            macros_file.write(self.cmake_string)
        with open(cmakelists_name, "w") as cmakelists:
            cmakelists.write(self._cmakelists_template.format("CIME_"+var_name))

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
        maker = Compilers(MockMachines("mymachine", "SomeOS"), version="2.0")
        test_xml = _wrap_config_build_xml("<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>")
        get_macros(maker, test_xml, "Makefile")

    def test_script_rejects_bad_xml(self):
        """The macro writer rejects input that's not valid XML."""
        maker = Compilers(MockMachines("mymachine", "SomeOS"), version="2.0")
        with self.assertRaises(ParseError):
            get_macros(maker, "This is not valid XML.", "Makefile")

    def test_script_rejects_bad_build_system(self):
        """The macro writer rejects a bad build system string."""
        maker = Compilers(MockMachines("mymachine", "SomeOS"), version="2.0")
        bad_string = "argle-bargle."
        with self.assertRaisesRegexp(
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

    test_os = "SomeOS"
    test_machine = "mymachine"

    def setUp(self):
        self._maker = Compilers(MockMachines(self.test_machine, self.test_os), version="2.0")

    def xml_to_tester(self, xml_string):
        """Helper that directly converts an XML string to a MakefileTester."""
        test_xml = _wrap_config_build_xml(xml_string)
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
        with self.assertRaisesRegexp(
                SystemExit,
                "Variable MPI_PATH is set ambiguously in config_build.xml."):
            self.xml_to_tester(xml1+xml2)

    def test_reject_duplicates(self):
        """The macro writer dies if given many matches for a given configuration."""
        xml1 = """<compiler><MPI_PATH MPILIB="mpich">/path/to/mpich</MPI_PATH></compiler>"""
        xml2 = """<compiler><MPI_PATH MPILIB="mpich">/path/to/mpich2</MPI_PATH></compiler>"""
        with self.assertRaisesRegexp(
                SystemExit,
                "Variable MPI_PATH is set ambiguously in config_build.xml."):
            self.xml_to_tester(xml1+xml2)

    def test_reject_ambiguous(self):
        """The macro writer dies if given an ambiguous set of matches."""
        xml1 = """<compiler><MPI_PATH MPILIB="mpich">/path/to/mpich</MPI_PATH></compiler>"""
        xml2 = """<compiler><MPI_PATH DEBUG="FALSE">/path/to/mpi-debug</MPI_PATH></compiler>"""
        with self.assertRaisesRegexp(
                SystemExit,
                "Variable MPI_PATH is set ambiguously in config_build.xml."):
            self.xml_to_tester(xml1+xml2)

    def test_compiler_changeable_at_build_time(self):
        """The macro writer writes information for multiple compilers."""
        xml1 = """<compiler><SUPPORTS_CXX>FALSE</SUPPORTS_CXX></compiler>"""
        xml2 = """<compiler COMPILER="gnu"><SUPPORTS_CXX>TRUE</SUPPORTS_CXX></compiler>"""
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("SUPPORTS_CXX", "FALSE")
        tester.assert_variable_equals("SUPPORTS_CXX", "TRUE", env={"COMPILER": "gnu"})

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

    def test_machine_specific_base_over_append_flags(self):
        """Test that machine-specific base flags override default append flags."""
        xml1 = """<compiler><FFLAGS><append>-delicious</append></FFLAGS></compiler>"""
        xml2 = """<compiler MACH="{}"><FFLAGS><base>-cake</base></FFLAGS></compiler>""".format(self.test_machine)
        tester = self.xml_to_tester(xml1+xml2)
        tester.assert_variable_equals("FFLAGS", "-cake")
        tester = self.xml_to_tester(xml2+xml1)
        tester.assert_variable_equals("FFLAGS", "-cake")

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
        """Test that <env> elements insert environment variables."""
        xml1 = """<compiler><LDFLAGS><append>-L<env>NETCDF</env> -lnetcdf</append></LDFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("LDFLAGS", "-L/path/to/netcdf -lnetcdf",
                                      env={"NETCDF": "/path/to/netcdf"})

    def test_shell_command_insertion(self):
        """Test that <shell> elements insert shell command output."""
        xml1 = """<compiler><FFLAGS><base>-O<shell>echo 2</shell> -fast</base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-O2 -fast")

    def test_multiple_shell_commands(self):
        """Test that more than one <shell> element can be used."""
        xml1 = """<compiler><FFLAGS><base>-O<shell>echo 2</shell> -<shell>echo fast</shell></base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-O2 -fast")

    def test_env_and_shell_command(self):
        """Test that <env> elements work inside <shell> elements."""
        xml1 = """<compiler><FFLAGS><base>-O<shell>echo <env>OPT_LEVEL</env></shell> -fast</base></FFLAGS></compiler>"""
        tester = self.xml_to_tester(xml1)
        tester.assert_variable_equals("FFLAGS", "-O2 -fast", env={"OPT_LEVEL": "2"})

    def test_config_variable_insertion(self):
        """Test that <var> elements insert variables from config_build."""
        # Construct an absurd chain of references just to sure that we don't
        # pass by accident, e.g. outputting things in the right order just due
        # to good luck in a hash somewhere.
        xml1 = """<MPI_LIB_NAME>stuff-<var>MPI_PATH</var>-stuff</MPI_LIB_NAME>"""
        xml2 = """<MPI_PATH><var>MPICC</var></MPI_PATH>"""
        xml3 = """<MPICC><var>MPICXX</var></MPICC>"""
        xml4 = """<MPICXX><var>MPIFC</var></MPICXX>"""
        xml5 = """<MPIFC>mpicc</MPIFC>"""
        tester = self.xml_to_tester("<compiler>"+xml1+xml2+xml3+xml4+xml5+"</compiler>")
        tester.assert_variable_equals("MPI_LIB_NAME", "stuff-mpicc-stuff")

    def test_config_reject_self_references(self):
        """Test that <var> self-references are rejected."""
        # This is a special case of the next test, which also checks circular
        # references.
        xml1 = """<MPI_LIB_NAME><var>MPI_LIB_NAME</var></MPI_LIB_NAME>"""
        err_msg = "The config_build XML has bad <var> references."
        with self.assertRaisesRegexp(SystemExit, err_msg):
            self.xml_to_tester("<compiler>"+xml1+"</compiler>")

    def test_config_reject_cyclical_references(self):
        """Test that cyclical <var> references are rejected."""
        xml1 = """<MPI_LIB_NAME><var>MPI_PATH</var></MPI_LIB_NAME>"""
        xml2 = """<MPI_PATH><var>MPI_LIB_NAME</var></MPI_PATH>"""
        err_msg = "The config_build XML has bad <var> references."
        with self.assertRaisesRegexp(SystemExit, err_msg):
            self.xml_to_tester("<compiler>"+xml1+xml2+"</compiler>")

    def test_variable_insertion_with_machine_specific_setting(self):
        """Test that machine-specific <var> dependencies are correct."""
        xml1 = """<compiler><MPI_LIB_NAME>something</MPI_LIB_NAME></compiler>"""
        xml2 = """<compiler MACH="{}"><MPI_LIB_NAME><var>MPI_PATH</var></MPI_LIB_NAME></compiler>""".format(self.test_machine)
        xml3 = """<compiler><MPI_PATH><var>MPI_LIB_NAME</var></MPI_PATH></compiler>"""
        err_msg = "The config_build XML has bad <var> references."
        with self.assertRaisesRegexp(SystemExit, err_msg):
            self.xml_to_tester(xml1+xml2+xml3)


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
        test_xml = _wrap_config_build_xml(xml_string)
        if (NO_CMAKE):
            self.skipTest("Skipping cmake test")
        else:
            return CMakeTester(self, get_macros(self._maker, test_xml, "CMake"))

###############################################################################
class S_TestManageAndQuery(unittest.TestCase):
    """Tests various scripts to manage and query xml files"""

    def _run_and_assert_query_testlist(self, extra_args=""):
        """Ensure that query_testlist runs successfully with the given extra arguments"""

        testlist_allactive = os.path.join(CIME.utils.get_model_config_root(), "allactive", "testlist_allactive.xml")

        run_cmd_assert_result(self, "%s/query_testlists --xml-testlist %s %s"%
                              (SCRIPT_DIR, testlist_allactive, extra_args))

    def test_query_testlists_runs(self):
        """Make sure that query_testlists runs successfully

        This simply makes sure that query_testlists doesn't generate any errors
        when it runs. This helps ensure that changes in other utilities don't
        break query_testlists.
        """
        self._run_and_assert_query_testlist(extra_args="--show-options")

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
        result = B_CheckCode.all_results[pyfile]
        self.assertTrue(result == "", msg=result)

    return test

def check_for_pylint():
    from distutils.spawn import find_executable
    pylint = find_executable("pylint")
    if pylint is not None:
        output = run_cmd_no_fail("pylint --version")
        pylintver = re.search(r"pylint\s+(\d+)[.](\d+)[.](\d+)", output)
        major = int(pylintver.group(1))
        minor = int(pylintver.group(2))
    if pylint is None or major < 1 or (major == 1 and minor < 5):
        print "pylint version 1.5 or newer not found, pylint tests skipped"
        return False
    return True

def write_provenance_info():
    curr_commit = get_current_commit(repo=LIB_DIR)
    logging.info("\nTesting commit %s" % curr_commit)
    cime_model = CIME.utils.get_model()
    logging.info("Using cime_model = %s\n" % cime_model)

def _main_func():

    if "--fast" in sys.argv:
        sys.argv.remove("--fast")
        global FAST_ONLY
        FAST_ONLY = True

    if "--no-batch" in sys.argv:
        sys.argv.remove("--no-batch")
        global NO_BATCH
        NO_BATCH = True

    if "--no-cmake" in sys.argv:
        sys.argv.remove("--no-cmake")
        global NO_CMAKE
        NO_CMAKE = True

    if "--test-root" in sys.argv:
        global TEST_ROOT
        trindex = sys.argv.index("--test-root")
        TEST_ROOT = sys.argv[trindex +1]
        del sys.argv[trindex+1]
        del sys.argv[trindex]
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

    CIME.utils.handle_standard_logging_options(args)

    write_provenance_info()

    # Find all python files in repo and create a pylint test for each
    if check_for_pylint():
        files_to_test = get_all_checkable_files()

        for file_to_test in files_to_test:
            pylint_test = make_pylint_test(file_to_test, files_to_test)
            testname = "test_pylint_%s" % file_to_test.replace("/", "_").replace(".", "_")
            expect(not hasattr(B_CheckCode, testname), "Repeat %s" % testname)
            setattr(B_CheckCode, testname, pylint_test)

    unittest.main(verbosity=2, catchbreak=True)

if (__name__ == "__main__"):
    _main_func()
