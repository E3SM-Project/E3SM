"""
A library for scheduling/running through the phases of a set
of system tests. Supports phase-level parallelism (can make progres
on multiple system tests at once).

TestScheduler will handle the TestStatus for the 1-time setup
phases. All other phases need to handle their own status because
they can be run outside the context of TestScheduler.
"""

import traceback, stat, threading, time, glob
from collections import OrderedDict

from CIME.XML.standard_module_setup import *
import six
from get_tests import get_recommended_test_time
from CIME.utils import append_status, append_testlog, TESTS_FAILED_ERR_CODE, parse_test_name, get_full_test_name, get_model, \
    convert_to_seconds, get_cime_root, get_project, get_timestamp, get_python_libs_root
from CIME.test_status import *
from CIME.XML.machines import Machines
from CIME.XML.generic_xml import GenericXML
from CIME.XML.env_test import EnvTest
from CIME.XML.env_mach_pes import EnvMachPes
from CIME.XML.files import Files
from CIME.XML.component import Component
from CIME.XML.tests import Tests
from CIME.case import Case
from CIME.wait_for_tests import wait_for_tests
from CIME.provenance import get_recommended_test_time_based_on_past
from CIME.locked_files import lock_file
from CIME.cs_status_creator import create_cs_status
from CIME.hist_utils import generate_teststatus

logger = logging.getLogger(__name__)

# Phases managed by TestScheduler
TEST_START = "INIT" # Special pseudo-phase just for test_scheduler bookkeeping
PHASES = [TEST_START, CREATE_NEWCASE_PHASE, XML_PHASE, SETUP_PHASE,
          SHAREDLIB_BUILD_PHASE, MODEL_BUILD_PHASE, RUN_PHASE] # Order matters

###############################################################################
def _translate_test_names_for_new_pecount(test_names, force_procs, force_threads):
###############################################################################
    new_test_names = []
    caseopts = []
    for test_name in test_names:
        testcase, caseopts, grid, compset, machine, compiler, testmod = parse_test_name(test_name)
        rewrote_caseopt = False
        if caseopts is not None:
            for idx, caseopt in enumerate(caseopts):
                if caseopt.startswith("P"):
                    caseopt = caseopt[1:]
                    if "x" in caseopt:
                        old_procs, old_thrds = caseopt.split("x")
                    else:
                        old_procs, old_thrds = caseopt, None

                    new_procs = force_procs if force_procs is not None else old_procs
                    new_thrds = force_threads if force_threads is not None else old_thrds

                    newcaseopt = ("P{}".format(new_procs)) if new_thrds is None else ("P{}x{}".format(new_procs, new_thrds))
                    caseopts[idx] = newcaseopt

                    rewrote_caseopt = True
                    break

        if not rewrote_caseopt:
            force_procs = "M" if force_procs is None else force_procs
            newcaseopt = ("P{}".format(force_procs)) if force_threads is None else ("P{}x{}".format(force_procs, force_threads))
            if caseopts is None:
                caseopts = [newcaseopt]
            else:
                caseopts.append(newcaseopt)

        new_test_name = get_full_test_name(testcase, caseopts=caseopts, grid=grid, compset=compset, machine=machine, compiler=compiler, testmod=testmod)
        new_test_names.append(new_test_name)

    return new_test_names

_TIME_CACHE = {}
###############################################################################
def _get_time_est(test, baseline_root, as_int=False, use_cache=False, raw=False):
###############################################################################
    if test in _TIME_CACHE and use_cache:
        return _TIME_CACHE[test]

    recommended_time = get_recommended_test_time_based_on_past(baseline_root, test, raw=raw)

    if recommended_time is None:
        recommended_time = get_recommended_test_time(test)

    if as_int:
        if recommended_time is None:
            recommended_time = 9999999999
        else:
            recommended_time = convert_to_seconds(recommended_time)

    if use_cache:
        _TIME_CACHE[test] = recommended_time

    return recommended_time

###############################################################################
def _order_tests_by_runtime(tests, baseline_root):
###############################################################################
    tests.sort(key=lambda x: _get_time_est(x, baseline_root, as_int=True, use_cache=True, raw=True), reverse=True)

###############################################################################
class TestScheduler(object):
###############################################################################

    ###########################################################################
    def __init__(self, test_names, test_data=None,
                 no_run=False, no_build=False, no_setup=False, no_batch=None,
                 test_root=None, test_id=None,
                 machine_name=None, compiler=None,
                 baseline_root=None, baseline_cmp_name=None, baseline_gen_name=None,
                 clean=False, namelists_only=False,
                 project=None, parallel_jobs=None,
                 walltime=None, proc_pool=None,
                 use_existing=False, save_timing=False, queue=None,
                 allow_baseline_overwrite=False, output_root=None,
                 force_procs=None, force_threads=None, mpilib=None,
                 input_dir=None, pesfile=None, mail_user=None, mail_type=None, allow_pnl=False, non_local=False):
    ###########################################################################
        self._cime_root       = get_cime_root()
        self._cime_model      = get_model()
        self._cime_driver     = "mct"
        self._save_timing     = save_timing
        self._queue           = queue
        self._test_data       = {} if test_data is None else test_data # Format:  {test_name -> {data_name -> data}}
        self._mpilib          = mpilib  # allow override of default mpilib
        self._completed_tests = 0
        self._input_dir       = input_dir
        self._pesfile         = pesfile
        self._allow_baseline_overwrite = allow_baseline_overwrite
        self._allow_pnl       = allow_pnl
        self._non_local       = non_local

        self._mail_user = mail_user
        self._mail_type = mail_type

        self._machobj = Machines(machine=machine_name)

        self._model_build_cost = 4

        # If user is forcing procs or threads, re-write test names to reflect this.
        if force_procs or force_threads:
            test_names = _translate_test_names_for_new_pecount(test_names, force_procs, force_threads)

        self._no_setup = no_setup
        self._no_build = no_build or no_setup or namelists_only
        self._no_run   = no_run or self._no_build
        self._output_root = output_root
        # Figure out what project to use
        if project is None:
            self._project = get_project()
            if self._project is None:
                self._project = self._machobj.get_value("PROJECT")
        else:
            self._project = project

        # We will not use batch system if user asked for no_batch or if current
        # machine is not a batch machine
        self._no_batch = no_batch or not self._machobj.has_batch_system()
        expect(not (self._no_batch and self._queue is not None),
               "Does not make sense to request a queue without batch system")

        # Determine and resolve test_root
        if test_root is not None:
            self._test_root = test_root
        elif self._output_root is not None:
            self._test_root = self._output_root
        else:
            self._test_root = self._machobj.get_value("CIME_OUTPUT_ROOT")

        if self._project is not None:
            self._test_root = self._test_root.replace("$PROJECT", self._project)

        self._test_root = os.path.abspath(self._test_root)
        self._test_id   = test_id if test_id is not None else get_timestamp()

        self._compiler = self._machobj.get_default_compiler() if compiler is None else compiler

        self._clean          = clean
        self._namelists_only = namelists_only

        self._walltime = walltime

        if parallel_jobs is None:
            self._parallel_jobs = min(len(test_names),
                                      self._machobj.get_value("MAX_MPITASKS_PER_NODE"))
        else:
            self._parallel_jobs = parallel_jobs

        self._baseline_cmp_name = baseline_cmp_name # Implies comparison should be done if not None
        self._baseline_gen_name = baseline_gen_name # Implies generation should be done if not None

        # Compute baseline_root
        self._baseline_root = baseline_root if baseline_root is not None \
                              else self._machobj.get_value("BASELINE_ROOT")

        if self._project is not None:
            self._baseline_root = self._baseline_root.replace("$PROJECT", self._project)

        self._baseline_root = os.path.abspath(self._baseline_root)

        if baseline_cmp_name or baseline_gen_name:
            if self._baseline_cmp_name:
                full_baseline_dir = os.path.join(self._baseline_root, self._baseline_cmp_name)
                expect(os.path.isdir(full_baseline_dir),
                       "Missing baseline comparison directory {}".format(full_baseline_dir))

            # the following is to assure that the existing generate directory is not overwritten
            if self._baseline_gen_name:
                full_baseline_dir = os.path.join(self._baseline_root, self._baseline_gen_name)
                existing_baselines = []
                for test_name in test_names:
                    test_baseline = os.path.join(full_baseline_dir, test_name)
                    if os.path.isdir(test_baseline):
                        existing_baselines.append(test_baseline)

                expect(allow_baseline_overwrite or len(existing_baselines) == 0,
                       "Baseline directories already exists {}\n" \
                       "Use -o to avoid this error".format(existing_baselines))

        if self._cime_model == "e3sm":
            _order_tests_by_runtime(test_names, self._baseline_root)

        # This is the only data that multiple threads will simultaneously access
        # Each test has it's own value and setting/retrieving items from a dict
        # is atomic, so this should be fine to use without mutex.
        # name -> (phase, status)
        self._tests = OrderedDict()
        for test_name in test_names:
            self._tests[test_name] = (TEST_START, TEST_PASS_STATUS)

        # Oversubscribe by 1/4
        if proc_pool is None:
            pes = int(self._machobj.get_value("MAX_TASKS_PER_NODE"))
            self._proc_pool = int(pes * 1.25)
        else:
            self._proc_pool = int(proc_pool)

        self._procs_avail = self._proc_pool

        # Setup phases
        self._phases = list(PHASES)
        if self._no_setup:
            self._phases.remove(SETUP_PHASE)
        if self._no_build:
            self._phases.remove(SHAREDLIB_BUILD_PHASE)
            self._phases.remove(MODEL_BUILD_PHASE)
        if self._no_run:
            self._phases.remove(RUN_PHASE)

        if use_existing:
            for test in self._tests:
                with TestStatus(self._get_test_dir(test)) as ts:
                    for phase, status in ts:
                        if phase in CORE_PHASES:
                            if status in [TEST_PEND_STATUS, TEST_FAIL_STATUS]:
                                if status == TEST_FAIL_STATUS:
                                    # Import for potential subsequent waits
                                    ts.set_status(phase, TEST_PEND_STATUS)

                                # We need to pick up here
                                break

                            else:
                                if phase != SUBMIT_PHASE:
                                    # Somewhat subtle. Create_test considers submit/run to be the run phase,
                                    # so don't try to update test status for a passed submit phase
                                    self._update_test_status(test, phase, TEST_PEND_STATUS)
                                    self._update_test_status(test, phase, status)

                                    if phase == RUN_PHASE:
                                        logger.info("Test {} passed and will not be re-run".format(test))

                logger.info("Using existing test directory {}".format(self._get_test_dir(test)))
        else:
            # None of the test directories should already exist.
            for test in self._tests:
                expect(not os.path.exists(self._get_test_dir(test)),
                       "Cannot create new case in directory '{}', it already exists."
                       " Pick a different test-id".format(self._get_test_dir(test)))
                logger.info("Creating test directory {}".format(self._get_test_dir(test)))

        # By the end of this constructor, this program should never hard abort,
        # instead, errors will be placed in the TestStatus files for the various
        # tests cases

    ###########################################################################
    def get_testnames(self):
    ###########################################################################
        return list(self._tests.keys())

    ###########################################################################
    def _log_output(self, test, output):
    ###########################################################################
        test_dir = self._get_test_dir(test)
        if not os.path.isdir(test_dir):
            # Note: making this directory could cause create_newcase to fail
            # if this is run before.
            os.makedirs(test_dir)
        append_testlog(output, caseroot=test_dir)

    ###########################################################################
    def _get_case_id(self, test):
    ###########################################################################
        baseline_action_code = ""
        if self._baseline_gen_name:
            baseline_action_code += "G"
        if self._baseline_cmp_name:
            baseline_action_code += "C"
        if len(baseline_action_code) > 0:
            return "{}.{}.{}".format(test, baseline_action_code, self._test_id)
        else:
            return "{}.{}".format(test, self._test_id)

    ###########################################################################
    def _get_test_dir(self, test):
    ###########################################################################
        return os.path.join(self._test_root, self._get_case_id(test))

    ###########################################################################
    def _get_test_data(self, test):
    ###########################################################################
        # Must be atomic
        return self._tests[test]

    ###########################################################################
    def _is_broken(self, test):
    ###########################################################################
        status = self._get_test_status(test)
        return status != TEST_PASS_STATUS and status != TEST_PEND_STATUS

    ###########################################################################
    def _work_remains(self, test):
    ###########################################################################
        test_phase, test_status = self._get_test_data(test)
        return (test_status == TEST_PASS_STATUS or test_status == TEST_PEND_STATUS) and\
            test_phase != self._phases[-1]

    ###########################################################################
    def _get_test_status(self, test, phase=None):
    ###########################################################################
        curr_phase, curr_status = self._get_test_data(test)
        if phase is None or phase == curr_phase:
            return curr_status
        else:
            expect(phase is None or self._phases.index(phase) < self._phases.index(curr_phase),
                   "Tried to see the future")
            # Assume all older phases PASSed
            return TEST_PASS_STATUS

    ###########################################################################
    def _get_test_phase(self, test):
    ###########################################################################
        return self._get_test_data(test)[0]

    ###########################################################################
    def _update_test_status(self, test, phase, status):
    ###########################################################################
        phase_idx = self._phases.index(phase)
        old_phase, old_status = self._get_test_data(test)

        if old_phase == phase:
            expect(old_status == TEST_PEND_STATUS,
                   "Only valid to transition from PEND to something else, found '{}' for phase '{}'".format(old_status, phase))
            expect(status != TEST_PEND_STATUS,
                   "Cannot transition from PEND -> PEND")
        else:
            expect(old_status == TEST_PASS_STATUS,
                   "Why did we move on to next phase when prior phase did not pass?")
            expect(status == TEST_PEND_STATUS,
                   "New phase should be set to pending status")
            expect(self._phases.index(old_phase) == phase_idx - 1,
                   "Skipped phase? {} {}".format(old_phase, phase_idx))

        # Must be atomic
        self._tests[test] = (phase, status)

    ###########################################################################
    def _shell_cmd_for_phase(self, test, cmd, phase, from_dir=None):
    ###########################################################################
        while True:
            rc, output, errput = run_cmd(cmd, from_dir=from_dir)
            if rc != 0:
                self._log_output(test,
                                 "{} FAILED for test '{}'.\nCommand: {}\nOutput: {}\n".
                                 format(phase, test, cmd,
                                        output.encode('utf-8') + b"\n" + errput.encode('utf-8')))
                # Temporary hack to get around odd file descriptor use by
                # buildnml scripts.
                if "bad interpreter" in output:
                    time.sleep(1)
                    continue
                else:
                    return False, errput
            else:
                # We don't want "RUN PASSED" in the TestStatus.log if the only thing that
                # succeeded was the submission.
                phase = "SUBMIT" if phase == RUN_PHASE else phase
                self._log_output(test,
                                 "{} PASSED for test '{}'.\nCommand: {}\nOutput: {}\n".
                                 format(phase, test, cmd,
                                        output.encode('utf-8') + b"\n" + errput.encode('utf-8')))
                return True, errput

    ###########################################################################
    def _create_newcase_phase(self, test):
    ###########################################################################
        test_dir = self._get_test_dir(test)

        _, case_opts, grid, compset,\
            machine, compiler, test_mods = parse_test_name(test)

        create_newcase_cmd = "{} --case {} --res {} --compset {}"\
                             " --test".format(os.path.join(self._cime_root, "scripts", "create_newcase"),
                                              test_dir, grid, compset)
        if machine is not None:
            create_newcase_cmd += " --machine {}".format(machine)
        if compiler is not None:
            create_newcase_cmd += " --compiler {}".format(compiler)
        if self._project is not None:
            create_newcase_cmd += " --project {} ".format(self._project)
        if self._output_root is not None:
            create_newcase_cmd += " --output-root {} ".format(self._output_root)
        if self._input_dir is not None:
            create_newcase_cmd += " --input-dir {} ".format(self._input_dir)
        if self._non_local:
            create_newcase_cmd += " --non-local"

        if self._pesfile is not None:
            create_newcase_cmd += " --pesfile {} ".format(self._pesfile)

        if test_mods is not None:
            files = Files(comp_interface=self._cime_driver)

            if test_mods.find('/') != -1:
                (component, modspath) = test_mods.split('/', 1)
            else:
                error = "Missing testmod component. Testmods are specified as '${component}-${testmod}'"
                self._log_output(test, error)
                return False, error

            testmods_dir = files.get_value("TESTS_MODS_DIR", {"component": component})
            test_mod_file = os.path.join(testmods_dir, component, modspath)
            if not os.path.exists(test_mod_file):
                error = "Missing testmod file '{}'".format(test_mod_file)
                self._log_output(test, error)
                return False, error

            create_newcase_cmd += " --user-mods-dir {}".format(test_mod_file)

        mpilib = None
        ninst = 1
        ncpl = 1
        if case_opts is not None:
            for case_opt in case_opts: # pylint: disable=not-an-iterable
                if case_opt.startswith('M'):
                    mpilib = case_opt[1:]
                    create_newcase_cmd += " --mpilib {}".format(mpilib)
                    logger.debug (" MPILIB set to {}".format(mpilib))
                elif case_opt.startswith('N'):
                    expect(ncpl == 1,"Cannot combine _C and _N options")
                    ninst = case_opt[1:]
                    create_newcase_cmd += " --ninst {}".format(ninst)
                    logger.debug (" NINST set to {}".format(ninst))
                elif case_opt.startswith('C'):
                    expect(ninst == 1,"Cannot combine _C and _N options")
                    ncpl = case_opt[1:]
                    create_newcase_cmd += " --ninst {} --multi-driver" .format(ncpl)
                    logger.debug (" NCPL set to {}" .format(ncpl))
                elif case_opt.startswith('P'):
                    pesize = case_opt[1:]
                    create_newcase_cmd += " --pecount {}".format(pesize)
                elif case_opt.startswith('V'):
                    self._cime_driver = case_opt[1:]
                    create_newcase_cmd += " --driver {}".format(self._cime_driver)


        # create_test mpilib option overrides default but not explicitly set case_opt mpilib
        if mpilib is None and self._mpilib is not None:
            create_newcase_cmd += " --mpilib {}".format(self._mpilib)
            logger.debug (" MPILIB set to {}".format(self._mpilib))

        if self._queue is not None:
            create_newcase_cmd += " --queue={}".format(self._queue)
        else:
            # We need to hard code the queue for this test on cheyenne
            # otherwise it runs in share and fails intermittently
            test_case = parse_test_name(test)[0]
            if test_case == "NODEFAIL":
                machine = machine if machine is not None else self._machobj.get_machine_name()
                if machine == "cheyenne":
                    create_newcase_cmd += " --queue=regular"

        if self._walltime is not None:
            create_newcase_cmd += " --walltime {}".format(self._walltime)
        else:
            # model specific ways of setting time
            if self._cime_model == "e3sm":
                recommended_time = _get_time_est(test, self._baseline_root)

                if recommended_time is not None:
                    create_newcase_cmd += " --walltime {}".format(recommended_time)

            else:
                if test in self._test_data and "options" in self._test_data[test] and \
                        "wallclock" in self._test_data[test]['options']:
                    create_newcase_cmd += " --walltime {}".format(self._test_data[test]['options']['wallclock'])

        logger.debug("Calling create_newcase: " + create_newcase_cmd)
        return self._shell_cmd_for_phase(test, create_newcase_cmd, CREATE_NEWCASE_PHASE)

    ###########################################################################
    def _xml_phase(self, test):
    ###########################################################################
        test_case = parse_test_name(test)[0]

        # Create, fill and write an envtest object
        test_dir = self._get_test_dir(test)
        envtest = EnvTest(test_dir)

        # Determine list of component classes that this coupler/driver knows how
        # to deal with. This list follows the same order as compset longnames follow.
        files = Files(comp_interface=self._cime_driver)
        drv_config_file = files.get_value("CONFIG_CPL_FILE")
        drv_comp = Component(drv_config_file, "CPL")
        envtest.add_elements_by_group(files, {}, "env_test.xml")
        envtest.add_elements_by_group(drv_comp, {}, "env_test.xml")
        envtest.set_value("TESTCASE", test_case)
        envtest.set_value("TEST_TESTID", self._test_id)
        envtest.set_value("CASEBASEID", test)
        if test in self._test_data and "options" in self._test_data[test] and \
                "memleak_tolerance" in self._test_data[test]['options']:
            envtest.set_value("TEST_MEMLEAK_TOLERANCE", self._test_data[test]['options']['memleak_tolerance'])

        test_argv = "-testname {} -testroot {}".format(test, self._test_root)
        if self._baseline_gen_name:
            test_argv += " -generate {}".format(self._baseline_gen_name)
            basegen_case_fullpath = os.path.join(self._baseline_root,self._baseline_gen_name, test)
            logger.debug("basegen_case is {}".format(basegen_case_fullpath))
            envtest.set_value("BASELINE_NAME_GEN", self._baseline_gen_name)
            envtest.set_value("BASEGEN_CASE", os.path.join(self._baseline_gen_name, test))
        if self._baseline_cmp_name:
            test_argv += " -compare {}".format(self._baseline_cmp_name)
            envtest.set_value("BASELINE_NAME_CMP", self._baseline_cmp_name)
            envtest.set_value("BASECMP_CASE", os.path.join(self._baseline_cmp_name, test))

        envtest.set_value("TEST_ARGV", test_argv)
        envtest.set_value("CLEANUP", self._clean)

        envtest.set_value("BASELINE_ROOT", self._baseline_root)
        envtest.set_value("GENERATE_BASELINE", self._baseline_gen_name is not None)
        envtest.set_value("COMPARE_BASELINE", self._baseline_cmp_name is not None)
        envtest.set_value("CCSM_CPRNC", self._machobj.get_value("CCSM_CPRNC", resolved=False))
        tput_tolerance = self._machobj.get_value("TEST_TPUT_TOLERANCE", resolved=False)
        if test in self._test_data and "options" in self._test_data[test] and \
           "tput_tolerance" in self._test_data[test]['options']:
            tput_tolerance = self._test_data[test]['options']['tput_tolerance']

        envtest.set_value("TEST_TPUT_TOLERANCE", 0.25 if tput_tolerance is None else tput_tolerance)

        # Add the test instructions from config_test to env_test in the case
        config_test = Tests()
        testnode = config_test.get_test_node(test_case)
        envtest.add_test(testnode)
        # Determine the test_case from the test name
        test_case, case_opts = parse_test_name(test)[:2]

        # Determine case_opts from the test_case
        if case_opts is not None:
            logger.debug("case_opts are {} ".format(case_opts))
            for opt in case_opts: # pylint: disable=not-an-iterable

                logger.debug("case_opt is {}".format(opt))
                if opt == 'D':
                    envtest.set_test_parameter("DEBUG", "TRUE")
                    logger.debug (" DEBUG set to TRUE")

                elif opt == 'E':
                    envtest.set_test_parameter("USE_ESMF_LIB", "TRUE")
                    logger.debug (" USE_ESMF_LIB set to TRUE")

                elif opt == 'CG':
                    envtest.set_test_parameter("CALENDAR", "GREGORIAN")
                    logger.debug (" CALENDAR set to {}".format(opt))

                elif opt.startswith('L'):
                    match =  re.match('L([A-Za-z])([0-9]*)', opt)
                    stop_option = {"y":"nyears", "m":"nmonths", "d":"ndays", "h":"nhours",
                                   "s":"nseconds", "n":"nsteps"}
                    opt = match.group(1)
                    envtest.set_test_parameter("STOP_OPTION",stop_option[opt])
                    opti = match.group(2)
                    envtest.set_test_parameter("STOP_N", opti)

                    logger.debug (" STOP_OPTION set to {}".format(stop_option[opt]))
                    logger.debug (" STOP_N      set to {}".format(opti))

                elif opt.startswith('R'):
                    # R option is for testing in PTS_MODE or Single Column Model
                    #  (SCM) mode
                    envtest.set_test_parameter("PTS_MODE", "TRUE")

                    # For PTS_MODE, compile with mpi-serial
                    envtest.set_test_parameter("MPILIB", "mpi-serial")

                elif (opt.startswith('I') or # Marker to distinguish tests with same name - ignored
                      opt.startswith('M') or # handled in create_newcase
                      opt.startswith('P') or # handled in create_newcase
                      opt.startswith('N') or # handled in create_newcase
                      opt.startswith('C') or # handled in create_newcase
                      opt.startswith('V')):  # handled in create_newcase
                    pass

                elif opt.startswith('IOP'):
                    logger.warning("IOP test option not yet implemented")
                else:
                    expect(False, "Could not parse option '{}' ".format(opt))

        envtest.write()
        lock_file("env_run.xml", caseroot=test_dir, newname="env_run.orig.xml")

        with Case(test_dir, read_only=False) as case:
            if self._output_root is None:
                self._output_root = case.get_value("CIME_OUTPUT_ROOT")
            # if we are running a single test we don't need sharedlibroot
            if len(self._tests) > 1 and self._cime_model != "e3sm":
                case.set_value("SHAREDLIBROOT",
                               os.path.join(self._output_root,
                                            "sharedlibroot.{}".format(self._test_id)))
            envtest.set_initial_values(case)
            case.set_value("TEST", True)
            case.set_value("SAVE_TIMING", self._save_timing)

            # Scale back build parallelism on systems with few cores
            if self._model_build_cost > self._proc_pool:
                case.set_value("GMAKE_J", self._proc_pool)
                self._model_build_cost = self._proc_pool

        return True, ""

    ###########################################################################
    def _setup_phase(self, test):
    ###########################################################################
        test_dir  = self._get_test_dir(test)
        rv = self._shell_cmd_for_phase(test, "./case.setup", SETUP_PHASE, from_dir=test_dir)

        # It's OK for this command to fail with baseline diffs but not catastrophically
        if rv[0]:
            cmdstat, output, _ = run_cmd("./case.cmpgen_namelists", combine_output=True, from_dir=test_dir)
            expect(cmdstat in [0, TESTS_FAILED_ERR_CODE], "Fatal error in case.cmpgen_namelists: {}".format(output))

        return rv

    ###########################################################################
    def _sharedlib_build_phase(self, test):
    ###########################################################################
        test_dir = self._get_test_dir(test)
        return self._shell_cmd_for_phase(test, "./case.build --sharedlib-only", SHAREDLIB_BUILD_PHASE, from_dir=test_dir)

    ###########################################################################
    def _model_build_phase(self, test):
    ###########################################################################
        test_dir = self._get_test_dir(test)
        return self._shell_cmd_for_phase(test, "./case.build --model-only", MODEL_BUILD_PHASE, from_dir=test_dir)

    ###########################################################################
    def _run_phase(self, test):
    ###########################################################################
        test_dir = self._get_test_dir(test)

        cmd = "./case.submit"
        if not self._allow_pnl:
            cmd += " --skip-preview-namelist"
        if self._no_batch:
            cmd += " --no-batch"
        if self._mail_user:
            cmd += " --mail-user={}".format(self._mail_user)
        if self._mail_type:
            cmd += " -M={}".format(",".join(self._mail_type))

        return self._shell_cmd_for_phase(test, cmd, RUN_PHASE, from_dir=test_dir)

    ###########################################################################
    def _run_catch_exceptions(self, test, phase, run):
    ###########################################################################
        try:
            return run(test)
        except (SystemExit, Exception) as e:
            exc_tb = sys.exc_info()[2]
            errput = "Test '{}' failed in phase '{}' with exception '{}'\n".format(test, phase, str(e))
            errput += ''.join(traceback.format_tb(exc_tb))
            self._log_output(test, errput)
            return False, errput

    ###########################################################################
    def _get_procs_needed(self, test, phase, threads_in_flight=None, no_batch=False):
    ###########################################################################
        if phase == RUN_PHASE and (self._no_batch or no_batch):
            test_dir = self._get_test_dir(test)
            total_pes = EnvMachPes(test_dir, read_only=True).get_value("TOTALPES")
            return total_pes

        elif (phase == SHAREDLIB_BUILD_PHASE):
            if self._cime_model == "cesm":
                # Will force serialization of sharedlib builds
                # TODO - instead of serializing, compute all library configs needed and build
                # them all in parallel
                for _, _, running_phase in threads_in_flight.values():
                    if (running_phase == SHAREDLIB_BUILD_PHASE):
                        return self._proc_pool + 1

            return 1
        elif (phase == MODEL_BUILD_PHASE):
            # Model builds now happen in parallel
            return self._model_build_cost
        else:
            return 1

    ###########################################################################
    def _wait_for_something_to_finish(self, threads_in_flight):
    ###########################################################################
        expect(len(threads_in_flight) <= self._parallel_jobs, "Oversubscribed?")
        finished_tests = []
        while not finished_tests:
            for test, thread_info in threads_in_flight.items():
                if not thread_info[0].is_alive():
                    finished_tests.append((test, thread_info[1]))

            if not finished_tests:
                time.sleep(0.2)

        for finished_test, procs_needed in finished_tests:
            self._procs_avail += procs_needed
            del threads_in_flight[finished_test]

    ###########################################################################
    def _update_test_status_file(self, test, test_phase, status):
    ###########################################################################
        """
        In general, test_scheduler should not be responsible for updating
        the TestStatus file, but there are a few cases where it has to.
        """
        test_dir = self._get_test_dir(test)
        with TestStatus(test_dir=test_dir, test_name=test) as ts:
            ts.set_status(test_phase, status)

    ###########################################################################
    def _consumer(self, test, test_phase, phase_method):
    ###########################################################################
        before_time = time.time()
        success, errors = self._run_catch_exceptions(test, test_phase, phase_method)
        elapsed_time = time.time() - before_time
        status  = (TEST_PEND_STATUS if test_phase == RUN_PHASE and not \
                   self._no_batch else TEST_PASS_STATUS) if success else TEST_FAIL_STATUS

        if status != TEST_PEND_STATUS:
            self._update_test_status(test, test_phase, status)

        if not self._work_remains(test):
            self._completed_tests += 1
            total = len(self._tests)
            status_str = "Finished {} for test {} in {:f} seconds ({}). [COMPLETED {:d} of {:d}]".format(test_phase, test, elapsed_time, status, self._completed_tests, total)
        else:
            status_str = "Finished {} for test {} in {:f} seconds ({})".format(test_phase, test, elapsed_time, status)

        if not success:
            status_str += "\n    Case dir: {}\n".format(self._get_test_dir(test))
            status_str += "    Errors were:\n        {}\n".format("\n        ".join(str(errors.encode('utf-8')).splitlines()))

        logger.info(status_str)

        if test_phase in [CREATE_NEWCASE_PHASE, XML_PHASE]:
            # These are the phases for which TestScheduler is reponsible for
            # updating the TestStatus file
            self._update_test_status_file(test, test_phase, status)

        if test_phase == XML_PHASE:
            append_status("Case Created using: "+" ".join(sys.argv), "README.case", caseroot=self._get_test_dir(test))

        # On batch systems, we want to immediately submit to the queue, because
        # it's very cheap to submit and will get us a better spot in line
        if (success and not self._no_run and not self._no_batch and test_phase == MODEL_BUILD_PHASE):
            logger.info("Starting {} for test {} with 1 proc on interactive node and {:d} procs on compute nodes".format(RUN_PHASE, test, self._get_procs_needed(test, RUN_PHASE, no_batch=True)))
            self._update_test_status(test, RUN_PHASE, TEST_PEND_STATUS)
            self._consumer(test, RUN_PHASE, self._run_phase)

    ###########################################################################
    def _producer(self):
    ###########################################################################
        threads_in_flight = {} # test-name -> (thread, procs, phase)
        while True:
            work_to_do = False
            num_threads_launched_this_iteration = 0
            for test in self._tests:
                logger.debug("test_name: " + test)

                if self._work_remains(test):
                    work_to_do = True

                    # If we have no workers available, immediately break out of loop so we can wait
                    if len(threads_in_flight) == self._parallel_jobs:
                        break

                    if test not in threads_in_flight:
                        test_phase, test_status = self._get_test_data(test)
                        expect(test_status != TEST_PEND_STATUS, test)
                        next_phase = self._phases[self._phases.index(test_phase) + 1]
                        procs_needed = self._get_procs_needed(test, next_phase, threads_in_flight)

                        if procs_needed <= self._procs_avail:
                            self._procs_avail -= procs_needed

                            # Necessary to print this way when multiple threads printing
                            logger.info("Starting {} for test {} with {:d} procs".format(next_phase, test, procs_needed))

                            self._update_test_status(test, next_phase, TEST_PEND_STATUS)
                            new_thread = threading.Thread(target=self._consumer,
                                args=(test, next_phase, getattr(self, "_{}_phase".format(next_phase.lower())) ))
                            threads_in_flight[test] = (new_thread, procs_needed, next_phase)
                            new_thread.start()
                            num_threads_launched_this_iteration += 1

                            logger.debug("  Current workload:")
                            total_procs = 0
                            for the_test, the_data in six.iteritems(threads_in_flight):
                                logger.debug("    {}: {} -> {}".format(the_test, the_data[2], the_data[1]))
                                total_procs += the_data[1]

                            logger.debug("    Total procs in use: {}".format(total_procs))
                        else:
                            if not threads_in_flight:
                                msg = "Phase '{}' for test '{}' required more processors, {:d}, than this machine can provide, {:d}".format(next_phase, test, procs_needed, self._procs_avail)
                                logger.warning(msg)
                                self._update_test_status(test, next_phase, TEST_PEND_STATUS)
                                self._update_test_status(test, next_phase, TEST_FAIL_STATUS)
                                self._log_output(test, msg)
                                if next_phase == RUN_PHASE:
                                    self._update_test_status_file(test, SUBMIT_PHASE, TEST_PASS_STATUS)
                                    self._update_test_status_file(test, next_phase, TEST_FAIL_STATUS)
                                else:
                                    self._update_test_status_file(test, next_phase, TEST_FAIL_STATUS)
                                num_threads_launched_this_iteration += 1

            if not work_to_do:
                break

            if num_threads_launched_this_iteration == 0:
                # No free resources, wait for something in flight to finish
                self._wait_for_something_to_finish(threads_in_flight)

        for unfinished_thread, _, _ in threads_in_flight.values():
            unfinished_thread.join()

    ###########################################################################
    def _setup_cs_files(self):
    ###########################################################################
        try:
            python_libs_root = get_python_libs_root()

            create_cs_status(test_root=self._test_root,
                             test_id=self._test_id)

            template_file = os.path.join(python_libs_root, "cs.submit.template")
            template = open(template_file, "r").read()
            setup_cmd = "./case.setup" if self._no_setup else ":"
            build_cmd = "./case.build" if self._no_build else ":"
            test_cmd  = "./case.submit"
            template = template.replace("<SETUP_CMD>", setup_cmd).\
                       replace("<BUILD_CMD>", build_cmd).\
                       replace("<RUN_CMD>", test_cmd).\
                       replace("<TESTID>", self._test_id)

            if self._no_run:
                cs_submit_file = os.path.join(self._test_root, "cs.submit.{}".format(self._test_id))
                with open(cs_submit_file, "w") as fd:
                    fd.write(template)
                os.chmod(cs_submit_file,
                         os.stat(cs_submit_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)

            if self._cime_model == "cesm":
                template_file = os.path.join(python_libs_root, "testreporter.template")
                template = open(template_file, "r").read()
                template = template.replace("<PATH>",
                                            os.path.join(self._cime_root, "scripts", "Tools"))
                testreporter_file = os.path.join(self._test_root, "testreporter")
                with open(testreporter_file, "w") as fd:
                    fd.write(template)
                os.chmod(testreporter_file, os.stat(testreporter_file).st_mode
                         | stat.S_IXUSR | stat.S_IXGRP)

        except Exception as e:
            logger.warning("FAILED to set up cs files: {}".format(str(e)))

    ###########################################################################
    def run_tests(self, wait=False,
                  wait_check_throughput=False,
                  wait_check_memory=False,
                  wait_ignore_namelists=False,
                  wait_ignore_memleak=False):
    ###########################################################################
        """
        Main API for this class.

        Return True if all tests passed.
        """
        start_time = time.time()

        # Tell user what will be run
        logger.info( "RUNNING TESTS:")
        for test in self._tests:
            logger.info( "  {}".format(test))

        # Setup cs files
        self._setup_cs_files()

        GenericXML.DISABLE_CACHING = True
        self._producer()
        GenericXML.DISABLE_CACHING = False

        expect(threading.active_count() == 1, "Leftover threads?")

        # Copy TestStatus files to baselines for tests that have already failed.
        if get_model() == "cesm":
            for test in self._tests:
                status = self._get_test_data(test)[1]
                if status not in [TEST_PASS_STATUS, TEST_PEND_STATUS] and self._baseline_gen_name:
                    basegen_case_fullpath = os.path.join(self._baseline_root,self._baseline_gen_name, test)
                    test_dir = self._get_test_dir(test)
                    generate_teststatus(test_dir, basegen_case_fullpath)

        wait_handles_report = False
        if not self._no_run and not self._no_batch:
            if wait:
                logger.info("Waiting for tests to finish")
                rv = wait_for_tests(glob.glob(os.path.join(self._test_root, "*{}/TestStatus".format(self._test_id))),
                                    check_throughput=wait_check_throughput,
                                    check_memory=wait_check_memory,
                                    ignore_namelists=wait_ignore_namelists,
                                    ignore_memleak=wait_ignore_memleak)
                wait_handles_report = True
            else:
                logger.info("Due to presence of batch system, create_test will exit before tests are complete.\n" \
                            "To force create_test to wait for full completion, use --wait")

        # Return True if all tests passed from our point of view
        if not wait_handles_report:
            logger.info( "At test-scheduler close, state is:")
            rv = True
            for test in self._tests:
                phase, status = self._get_test_data(test)

                # Give highest priority to fails in test schduler
                if status not in [TEST_PASS_STATUS, TEST_PEND_STATUS]:
                    logger.info( "{} {} (phase {})".format(status, test, phase))
                    rv = False

                else:
                    # Be cautious about telling the user that the test passed. This
                    # status should match what they would see on the dashboard. Our
                    # self._test_states does not include comparison fail information,
                    # so we need to parse test status.
                    ts = TestStatus(self._get_test_dir(test))
                    nlfail = ts.get_status(NAMELIST_PHASE) == TEST_FAIL_STATUS
                    ts_status = ts.get_overall_test_status(ignore_namelists=True, check_memory=False, check_throughput=False)
                    local_run = not self._no_run and self._no_batch

                    if ts_status not in [TEST_PASS_STATUS, TEST_PEND_STATUS]:
                        logger.info( "{} {} (phase {})".format(ts_status, test, phase))
                        rv = False
                    elif ts_status == TEST_PEND_STATUS and local_run:
                        logger.info( "{} {} (Some phases left in PEND)".format(TEST_FAIL_STATUS, test))
                        rv = False
                    elif nlfail:
                        logger.info( "{} {} (but otherwise OK) {}".format(NAMELIST_FAIL_STATUS, test, phase))
                        rv = False
                    else:
                        logger.info("{} {} {}".format(status, test, phase))

                logger.info( "    Case dir: {}".format(self._get_test_dir(test)))

        logger.info( "test-scheduler took {} seconds".format(time.time() - start_time))

        return rv
