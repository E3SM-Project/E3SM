"""
A library for scheduling/running through the phases of a set
of system tests. Supports phase-level parallelism (can make progres
on multiple system tests at once).

TestScheduler will handle the TestStatus for the 1-time setup
phases. All other phases need to handle their own status because
they can be run outside the context of TestScheduler.
"""

import shutil, traceback, stat, threading, time, glob
from CIME.XML.standard_module_setup import *
import CIME.compare_namelists
import CIME.utils
from CIME.utils import append_status, TESTS_FAILED_ERR_CODE, parse_test_name, get_full_test_name
from CIME.test_status import *
from CIME.XML.machines import Machines
from CIME.XML.env_test import EnvTest
from CIME.XML.files import Files
from CIME.XML.component import Component
from CIME.XML.tests import Tests
from CIME.case import Case
from CIME.wait_for_tests import wait_for_tests
from CIME.check_lockedfiles import lock_file
import CIME.test_utils

logger = logging.getLogger(__name__)

# Phases managed by TestScheduler
TEST_START = "INIT" # Special pseudo-phase just for test_scheduler bookkeeping
PHASES = [TEST_START, CREATE_NEWCASE_PHASE, XML_PHASE, SETUP_PHASE,
          SHAREDLIB_BUILD_PHASE, MODEL_BUILD_PHASE, RUN_PHASE] # Order matters

###############################################################################
def _translate_test_names_for_new_pecount(test_names, force_procs, force_threads):
###############################################################################
    new_test_names = []
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

                    newcaseopt = ("P%s" % new_procs) if new_thrds is None else ("P%sx%s" % (new_procs, new_thrds))

                    # No idea why pylint thinks this is unsubscriptable
                    caseopts[idx] = newcaseopt # pylint: disable=unsubscriptable-object

                    rewrote_caseopt = True
                    break

        if not rewrote_caseopt:
            force_procs = "M" if force_procs is None else force_procs
            newcaseopt = ("P%s" % force_procs) if force_threads is None else ("P%sx%s" % (force_procs, force_threads))
            if caseopts is None:
                caseopts = [newcaseopt]
            else:
                caseopts.append(newcaseopt)

        new_test_name = get_full_test_name(testcase, caseopts=caseopts, grid=grid, compset=compset, machine=machine, compiler=compiler, testmod=testmod)
        new_test_names.append(new_test_name)

    return new_test_names

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
                 force_procs=None, force_threads=None):
    ###########################################################################
        self._cime_root     = CIME.utils.get_cime_root()
        self._cime_model    = CIME.utils.get_model()
        self._save_timing   = save_timing
        self._queue         = queue
        self._test_data     = {} if test_data is None else test_data # Format:  {test_name -> {data_name -> data}}

        self._allow_baseline_overwrite  = allow_baseline_overwrite

        self._machobj = Machines(machine=machine_name)

        # If user is forcing procs or threads, re-write test names to reflect this.
        if force_procs or force_threads:
            test_names = _translate_test_names_for_new_pecount(test_names, force_procs, force_threads)

        self._no_setup = no_setup
        self._no_build = no_build or no_setup or namelists_only
        self._no_run   = no_run or self._no_build
        self._output_root = output_root
        # Figure out what project to use
        if project is None:
            self._project = CIME.utils.get_project()
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
        self._test_id   = test_id if test_id is not None else CIME.utils.get_timestamp()

        self._compiler = self._machobj.get_default_compiler() if compiler is None else compiler

        self._clean          = clean
        self._namelists_only = namelists_only

        self._walltime = walltime

        if parallel_jobs is None:
            self._parallel_jobs = min(len(test_names),
                                      int(self._machobj.get_value("MAX_TASKS_PER_NODE")))
        else:
            self._parallel_jobs = parallel_jobs

        self._baseline_cmp_name = baseline_cmp_name # Implies comparison should be done if not None
        self._baseline_gen_name = baseline_gen_name # Implies generation should be done if not None

        if baseline_cmp_name or baseline_gen_name:
            # Compute baseline_root
            self._baseline_root = baseline_root if baseline_root is not None \
                else self._machobj.get_value("BASELINE_ROOT")

            if self._project is not None:
                self._baseline_root = self._baseline_root.replace("$PROJECT", self._project)

            self._baseline_root = os.path.abspath(self._baseline_root)

            if self._baseline_cmp_name:
                full_baseline_dir = os.path.join(self._baseline_root, self._baseline_cmp_name)
                expect(os.path.isdir(full_baseline_dir),
                       "Missing baseline comparison directory %s" % full_baseline_dir)

            # the following is to assure that the existing generate directory is not overwritten
            if self._baseline_gen_name:
                full_baseline_dir = os.path.join(self._baseline_root, self._baseline_gen_name)
                existing_baselines = []
                for test_name in test_names:
                    test_baseline = os.path.join(full_baseline_dir, test_name)
                    if os.path.isdir(test_baseline):
                        existing_baselines.append(test_baseline)
                expect(allow_baseline_overwrite or len(existing_baselines) == 0,
                           "Baseline directories already exists %s\n"\
                           "Use --allow_baseline_overwrite to avoid this error"%existing_baselines)
        else:
            self._baseline_root = None

        # This is the only data that multiple threads will simultaneously access
        # Each test has it's own value and setting/retrieving items from a dict
        # is atomic, so this should be fine to use without mutex.
        # name -> (phase, status)
        self._tests = {}
        for test_name in test_names:
            self._tests[test_name] = (TEST_START, TEST_PASS_STATUS)

        # Oversubscribe by 1/4
        if proc_pool is None:
            pes = int(self._machobj.get_value("PES_PER_NODE"))
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
                ts = TestStatus(self._get_test_dir(test))
                for phase, status in ts:
                    if phase in CORE_PHASES:
                        if status in [TEST_PEND_STATUS, TEST_FAIL_STATUS]:
                            # We need to pick up here
                            break
                        else:
                            self._update_test_status(test, phase, TEST_PEND_STATUS)
                            self._update_test_status(test, phase, status)
        else:
            # None of the test directories should already exist.
            for test in self._tests:
                expect(not os.path.exists(self._get_test_dir(test)),
                       "Cannot create new case in directory '%s', it already exists."
                       " Pick a different test-id" % self._get_test_dir(test))
        logger.info("Created test in directory %s"%self._get_test_dir(test))
        # By the end of this constructor, this program should never hard abort,
        # instead, errors will be placed in the TestStatus files for the various
        # tests cases

    ###########################################################################
    def _log_output(self, test, output):
    ###########################################################################
        test_dir = self._get_test_dir(test)
        if not os.path.isdir(test_dir):
            # Note: making this directory could cause create_newcase to fail
            # if this is run before.
            os.makedirs(test_dir)
        append_status(output,caseroot=test_dir,sfile="TestStatus.log")

    ###########################################################################
    def _get_case_id(self, test):
    ###########################################################################
        baseline_action_code = ""
        if self._baseline_gen_name:
            baseline_action_code += "G"
        if self._baseline_cmp_name:
            baseline_action_code += "C"
        if len(baseline_action_code) > 0:
            return "%s.%s.%s" % (test, baseline_action_code, self._test_id)
        else:
            return "%s.%s" % (test, self._test_id)

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
                   "Only valid to transition from PEND to something else, found '%s' for phase '%s'" %
                   (old_status, phase))
            expect(status != TEST_PEND_STATUS,
                   "Cannot transition from PEND -> PEND")
        else:
            expect(old_status == TEST_PASS_STATUS,
                   "Why did we move on to next phase when prior phase did not pass?")
            expect(status == TEST_PEND_STATUS,
                   "New phase should be set to pending status")
            expect(self._phases.index(old_phase) == phase_idx - 1,
                   "Skipped phase? %s %s"%(old_phase, phase_idx))

        # Must be atomic
        self._tests[test] = (phase, status)

    ###########################################################################
    def _shell_cmd_for_phase(self, test, cmd, phase, from_dir=None):
    ###########################################################################
        while True:
            rc, output, errput = run_cmd(cmd, from_dir=from_dir)
            if rc != 0:
                self._log_output(test,
                                 "%s FAILED for test '%s'.\nCommand: %s\nOutput: %s\n\nErrput: %s" %
                                 (phase, test, cmd, output, errput))
                # Temporary hack to get around odd file descriptor use by
                # buildnml scripts.
                if "bad interpreter" in errput:
                    time.sleep(1)
                    continue
                else:
                    break
            else:
                # We don't want "RUN PASSED" in the TestStatus.log if the only thing that
                # succeeded was the submission.
                if phase != RUN_PHASE or self._no_batch:
                    self._log_output(test,
                                     "%s PASSED for test '%s'.\nCommand: %s\nOutput: %s" %
                                     (phase, test, cmd, output))
                break

        return rc == 0

    ###########################################################################
    def _create_newcase_phase(self, test):
    ###########################################################################
        test_dir = self._get_test_dir(test)

        _, case_opts, grid, compset,\
            machine, compiler, test_mods = CIME.utils.parse_test_name(test)

        create_newcase_cmd = "%s --case %s --res %s --mach %s --compiler %s --compset %s"\
                               " --test" % \
                              (os.path.join(self._cime_root, "scripts", "create_newcase"),
                               test_dir, grid, machine, compiler, compset)
        if self._project is not None:
            create_newcase_cmd += " --project %s " % self._project
        if self._output_root is not None:
            create_newcase_cmd += " --output-root %s " % self._output_root

        if test_mods is not None:
            files = Files()
            (component,modspath) = test_mods.split('/',1)
            testmods_dir = files.get_value("TESTS_MODS_DIR", {"component": component})

            test_mod_file = os.path.join(testmods_dir, component, modspath)
            if not os.path.exists(test_mod_file):
                self._log_output(test, "Missing testmod file '%s'" % test_mod_file)
                return False
            create_newcase_cmd += " --user-mods-dir %s" % test_mod_file

        if case_opts is not None:
            for case_opt in case_opts: # pylint: disable=not-an-iterable
                if case_opt.startswith('M'):
                    mpilib = case_opt[1:]
                    create_newcase_cmd += " --mpilib %s" % mpilib
                    logger.debug (" MPILIB set to %s" % mpilib)
                if case_opt.startswith('N'):
                    ninst = case_opt[1:]
                    create_newcase_cmd += " --ninst %s" %ninst
                    logger.debug (" NINST set to %s" % ninst)
                if case_opt.startswith('P'):
                    pesize = case_opt[1:]
                    create_newcase_cmd += " --pecount %s"%pesize

        if self._queue is not None:
            create_newcase_cmd += " --queue=%s" % self._queue

        if self._walltime is not None:
            create_newcase_cmd += " --walltime %s" % self._walltime
        elif test in self._test_data and "options" in self._test_data[test] and \
                "wallclock" in self._test_data[test]['options']:
            create_newcase_cmd += " --walltime %s" % self._test_data[test]['options']['wallclock']

        logger.debug("Calling create_newcase: " + create_newcase_cmd)
        return self._shell_cmd_for_phase(test, create_newcase_cmd, CREATE_NEWCASE_PHASE)

    ###########################################################################
    def _xml_phase(self, test):
    ###########################################################################
        test_case = CIME.utils.parse_test_name(test)[0]

        # Create, fill and write an envtest object
        test_dir = self._get_test_dir(test)
        envtest = EnvTest(test_dir)

        # Determine list of component classes that this coupler/driver knows how
        # to deal with. This list follows the same order as compset longnames follow.
        files = Files()
        drv_config_file = files.get_value("CONFIG_CPL_FILE")
        drv_comp = Component(drv_config_file)
        envtest.add_elements_by_group(files, {}, "env_test.xml")
        envtest.add_elements_by_group(drv_comp, {}, "env_test.xml")
        envtest.set_value("TESTCASE", test_case)
        envtest.set_value("TEST_TESTID", self._test_id)
        envtest.set_value("CASEBASEID", test)
        if test in self._test_data and "options" in self._test_data[test] and \
                "memleak_tolerance" in self._test_data[test]['options']:
            envtest.set_value("TEST_MEMLEAK_TOLERANCE", self._test_data[test]['options']['memleak_tolerance'])

        test_argv = "-testname %s -testroot %s" % (test, self._test_root)
        if self._baseline_gen_name:
            test_argv += " -generate %s" % self._baseline_gen_name
            basegen_case_fullpath = os.path.join(self._baseline_root,self._baseline_gen_name, test)
            logger.debug("basegen_case is %s"%basegen_case_fullpath)
            envtest.set_value("BASELINE_NAME_GEN", self._baseline_gen_name)
            envtest.set_value("BASEGEN_CASE", os.path.join(self._baseline_gen_name, test))
        if self._baseline_cmp_name:
            test_argv += " -compare %s" % self._baseline_cmp_name
            envtest.set_value("BASELINE_NAME_CMP", self._baseline_cmp_name)
            envtest.set_value("BASECMP_CASE", os.path.join(self._baseline_cmp_name, test))

        envtest.set_value("TEST_ARGV", test_argv)
        envtest.set_value("CLEANUP", self._clean)

        if self._baseline_gen_name or self._baseline_cmp_name:
            envtest.set_value("BASELINE_ROOT", self._baseline_root)
        envtest.set_value("GENERATE_BASELINE", self._baseline_gen_name is not None)
        envtest.set_value("COMPARE_BASELINE", self._baseline_cmp_name is not None)
        envtest.set_value("CCSM_CPRNC", self._machobj.get_value("CCSM_CPRNC", resolved=False))

        # Add the test instructions from config_test to env_test in the case
        config_test = Tests()
        testnode = config_test.get_test_node(test_case)
        envtest.add_test(testnode)

        # Determine the test_case from the test name
        test_case, case_opts = CIME.utils.parse_test_name(test)[:2]

        # Determine case_opts from the test_case
        if case_opts is not None:
            logger.debug("case_opts are %s " %case_opts)
            for opt in case_opts:

                logger.debug("case_opt is %s" %opt)
                if opt == 'D':
                    envtest.set_test_parameter("DEBUG", "TRUE")
                    logger.debug (" DEBUG set to TRUE")

                elif opt == 'E':
                    envtest.set_test_parameter("USE_ESMF_LIB", "TRUE")
                    envtest.set_test_parameter("COMP_INTERFACE", "ESMF")
                    logger.debug (" USE_ESMF_LIB set to TRUE")
                    logger.debug (" COMP_INTERFACE set to ESMF")

                elif opt == 'CG':
                    envtest.set_test_parameter("CALENDAR", "GREGORIAN")
                    logger.debug (" CALENDAR set to %s" %opt)

                elif opt.startswith('L'):
                    match =  re.match('L([A-Za-z])([0-9]*)', opt)
                    stop_option = {"y":"nyears", "m":"nmonths", "d":"ndays", "h":"nhours",
                                   "s":"nseconds", "n":"nsteps"}
                    opt = match.group(1)
                    envtest.set_test_parameter("STOP_OPTION",stop_option[opt])
                    opti = match.group(2)
                    envtest.set_test_parameter("STOP_N", opti)
                    logger.debug (" STOP_OPTION set to %s" %stop_option[opt])
                    logger.debug (" STOP_N      set to %s" %opti)

                elif opt.startswith('M'):
                    # M option handled by create newcase
                    continue

                elif opt.startswith('P'):
                    # P option handled by create newcase
                    continue

                elif opt.startswith('N'):
                    # handled in create_newcase
                    continue
                elif opt.startswith('IOP'):
                    logger.warn("IOP test option not yet implemented")
                else:
                    expect(False, "Could not parse option '%s' " %opt)

        envtest.write()
        lock_file("env_run.xml", caseroot=test_dir, newname="env_run.orig.xml")

        with Case(test_dir, read_only=False) as case:
            if self._output_root is None:
                self._output_root = case.get_value("CIME_OUTPUT_ROOT")
            # if we are running a single test we don't need sharedlibroot
            if len(self._tests) > 1:
                case.set_value("SHAREDLIBROOT",
                               os.path.join(self._output_root,
                                            "sharedlibroot.%s"%self._test_id))
            envtest.set_initial_values(case)
            case.set_value("TEST", True)
            if self._save_timing:
                case.set_value("SAVE_TIMING", True)

        return True

    ###########################################################################
    def _setup_phase(self, test):
    ###########################################################################
        test_dir  = self._get_test_dir(test)
        rv = self._shell_cmd_for_phase(test, "./case.setup", SETUP_PHASE, from_dir=test_dir)

        # A little subtle. If namelists_only, the RUN phase, when the namelists would normally
        # be handled, is not going to happen, so we have to do it here.
        if self._namelists_only:
            # It's OK for this command to fail with baseline diffs but not catastrophically
            cmdstat, output, errput = run_cmd("./case.cmpgen_namelists", from_dir=test_dir)
            expect(cmdstat in [0, TESTS_FAILED_ERR_CODE], "Fatal error in case.cmpgen_namelists: %s" % (output + "\n" + errput))

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
        if self._no_batch:
            cmd = "./case.submit --no-batch"
        else:
            cmd = "./case.submit "

        return self._shell_cmd_for_phase(test, cmd, RUN_PHASE, from_dir=test_dir)

    ###########################################################################
    def _run_catch_exceptions(self, test, phase, run):
    ###########################################################################
        try:
            return run(test)
        except (SystemExit, Exception) as e:
            exc_tb = sys.exc_info()[2]
            errput = "Test '%s' failed in phase '%s' with exception '%s'" % (test, phase, str(e))
            self._log_output(test, errput)
            logger.warning("Caught exception: %s" % str(e))
            traceback.print_tb(exc_tb)
            return False

    ###########################################################################
    def _get_procs_needed(self, test, phase, threads_in_flight=None, no_batch=False):
    ###########################################################################
        if phase == RUN_PHASE and (self._no_batch or no_batch):
            test_dir = self._get_test_dir(test)
            out = run_cmd_no_fail("./xmlquery TOTAL_CORES -value", from_dir=test_dir)
            return int(out)
        elif (phase == SHAREDLIB_BUILD_PHASE):
            # Will force serialization of sharedlib builds
            # TODO - instead of serializing, compute all library configs needed and build
            # them all in parallel
            for _, _, running_phase in threads_in_flight.values():
                if (running_phase == SHAREDLIB_BUILD_PHASE):
                    return self._proc_pool + 1

            return 1
        elif (phase == MODEL_BUILD_PHASE):
            # Model builds now happen in parallel
            return 4
        else:
            return 1

    ###########################################################################
    def _wait_for_something_to_finish(self, threads_in_flight):
    ###########################################################################
        expect(len(threads_in_flight) <= self._parallel_jobs, "Oversubscribed?")
        finished_tests = []
        while not finished_tests:
            for test, thread_info in threads_in_flight.iteritems():
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
        success = self._run_catch_exceptions(test, test_phase, phase_method)
        elapsed_time = time.time() - before_time
        status  = (TEST_PEND_STATUS if test_phase == RUN_PHASE and not \
                   self._no_batch else TEST_PASS_STATUS) if success else TEST_FAIL_STATUS

        if status != TEST_PEND_STATUS:
            self._update_test_status(test, test_phase, status)

        status_str = "Finished %s for test %s in %f seconds (%s)" %\
                     (test_phase, test, elapsed_time, status)
        if not success:
            status_str += "    Case dir: %s" % self._get_test_dir(test)
        logger.info(status_str)

        if test_phase in [CREATE_NEWCASE_PHASE, XML_PHASE]:
            # These are the phases for which TestScheduler is reponsible for
            # updating the TestStatus file
            self._update_test_status_file(test, test_phase, status)

        if test_phase == XML_PHASE:
            append_status("Case Created using: "+" ".join(sys.argv), caseroot=self._get_test_dir(test), sfile="README.case")

        # On batch systems, we want to immediately submit to the queue, because
        # it's very cheap to submit and will get us a better spot in line
        if (success and not self._no_run and not self._no_batch and test_phase == MODEL_BUILD_PHASE):
            logger.info("Starting %s for test %s with 1 proc on interactive node and %d procs on compute nodes" %
                        (RUN_PHASE, test, self._get_procs_needed(test, RUN_PHASE, no_batch=True)))
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
                # If we have no workers available, immediately wait
                if len(threads_in_flight) == self._parallel_jobs:
                    self._wait_for_something_to_finish(threads_in_flight)

                if self._work_remains(test):
                    work_to_do = True
                    if test not in threads_in_flight:
                        test_phase, test_status = self._get_test_data(test)
                        expect(test_status != TEST_PEND_STATUS, test)
                        next_phase = self._phases[self._phases.index(test_phase) + 1]
                        procs_needed = self._get_procs_needed(test, next_phase, threads_in_flight)

                        if procs_needed <= self._procs_avail:
                            self._procs_avail -= procs_needed

                            # Necessary to print this way when multiple threads printing
                            logger.info("Starting %s for test %s with %d procs" %
                                        (next_phase, test, procs_needed))

                            self._update_test_status(test, next_phase, TEST_PEND_STATUS)
                            new_thread = threading.Thread(target=self._consumer,
                                args=(test, next_phase, getattr(self, "_%s_phase" % next_phase.lower())) )
                            threads_in_flight[test] = (new_thread, procs_needed, next_phase)
                            new_thread.start()
                            num_threads_launched_this_iteration += 1
                        else:
                            if not threads_in_flight:
                                msg = "Phase '%s' for test '%s' required more processors, %d, than this machine can provide, %d" % \
                                    (next_phase, test, procs_needed, self._procs_avail)
                                logger.warning(msg)
                                self._update_test_status(test, next_phase, TEST_PEND_STATUS)
                                self._update_test_status(test, next_phase, TEST_FAIL_STATUS)
                                self._log_output(test, msg)
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
            python_libs_root = CIME.utils.get_python_libs_root()
            template_file = os.path.join(python_libs_root, "cs.status.template")
            template = open(template_file, "r").read()
            template = template.replace("<PATH>",
                                        os.path.join(self._cime_root,"scripts","Tools")).replace\
                                        ("<TESTID>", self._test_id)
            cs_status_file = os.path.join(self._test_root, "cs.status.%s" % self._test_id)
            with open(cs_status_file, "w") as fd:
                fd.write(template)
            os.chmod(cs_status_file, os.stat(cs_status_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)

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
                cs_submit_file = os.path.join(self._test_root, "cs.submit.%s" % self._test_id)
                with open(cs_submit_file, "w") as fd:
                    fd.write(template)
                os.chmod(cs_submit_file,
                         os.stat(cs_submit_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)

            if CIME.utils.get_model == "cesm":
                testreporter =  os.path.join(self._test_root,"testreporter.pl")
                shutil.copy(os.path.join(self._cime_root,"scripts","Testing","testreporter.pl"),
                            testreporter)
                os.chmod(testreporter, os.stat(testreporter).st_mode | stat.S_IXUSR | stat.S_IXGRP)

        except Exception as e:
            logger.warning("FAILED to set up cs files: %s" % str(e))

    ###########################################################################
    def run_tests(self, wait=False):
    ###########################################################################
        """
        Main API for this class.

        Return True if all tests passed.
        """
        start_time = time.time()

        # Tell user what will be run
        logger.info( "RUNNING TESTS:")
        for test in self._tests:
            logger.info( "  %s"% test)

        self._producer()

        expect(threading.active_count() == 1, "Leftover threads?")

        # Setup cs files
        self._setup_cs_files()

        wait_handles_report = False
        if not self._no_run and not self._no_batch:
            if wait:
                logger.info("Waiting for tests to finish")
                rv = wait_for_tests(glob.glob(os.path.join(self._test_root, "*%s/TestStatus" % self._test_id)))
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
                    logger.info( "%s %s (phase %s)" % (status, test, phase))
                    rv = False

                else:
                    # Be cautious about telling the user that the test passed. This
                    # status should match what they would see on the dashboard. Our
                    # self._test_states does not include comparison fail information,
                    # so we need to parse test status.
                    ts = TestStatus(self._get_test_dir(test))
                    nlfail = ts.get_status(NAMELIST_PHASE) == TEST_FAIL_STATUS
                    ts_status = ts.get_overall_test_status(ignore_namelists=True, check_memory=False, check_throughput=False)

                    if ts_status not in [TEST_PASS_STATUS, TEST_PEND_STATUS]:
                        logger.info( "%s %s (phase %s)" % (status, test, phase))
                        rv = False
                    elif nlfail:
                        logger.info( "%s %s (but otherwise OK) %s" % (NAMELIST_FAIL_STATUS, test, phase))
                        rv = False
                    else:
                        logger.info("%s %s %s" % (status, test, phase))

                logger.info( "    Case dir: %s" % self._get_test_dir(test))

            logger.info( "test-scheduler took %s seconds"% (time.time() - start_time))

        return rv
