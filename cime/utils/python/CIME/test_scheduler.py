"""
A library for scheduling/running through the phases of a set
of system tests. Supports phase-level parallelism (can make progres
on multiple system tests at once).

TestScheduler will handle the TestStatus for the 1-time setup
phases. All other phases need to handle their own status because
they can be run outside the context of TestScheduler.
"""

import shutil, traceback, stat, glob, threading, time
from CIME.XML.standard_module_setup import *
import CIME.compare_namelists
import CIME.utils
from CIME.utils import append_status
from CIME.test_status import *
from CIME.XML.machines import Machines
from CIME.XML.env_test import EnvTest
from CIME.XML.files import Files
from CIME.XML.component import Component
from CIME.XML.tests import Tests
from CIME.case import Case
import CIME.test_utils

logger = logging.getLogger(__name__)

# Phases managed by TestScheduler
PHASES = [INITIAL_PHASE, CREATE_NEWCASE_PHASE, XML_PHASE, SETUP_PHASE,
          NAMELIST_PHASE, SHAREDLIB_BUILD_PHASE, MODEL_BUILD_PHASE, RUN_PHASE] # Order matters
CONTINUE = [TEST_PASS_STATUS, NAMELIST_FAIL_STATUS]

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
                 use_existing=False, save_timing=False, queue=None, allow_baseline_overwrite=False):
    ###########################################################################
        self._cime_root  = CIME.utils.get_cime_root()
        self._cime_model = CIME.utils.get_model()
        self._allow_baseline_overwrite  = allow_baseline_overwrite
        self._save_timing = save_timing
        self._queue       = queue
        self._test_data   = {} if test_data is None else test_data # Format:  {test_name -> {data_name -> data}}

        # needed for perl interface
        os.environ["CIMEROOT"] = self._cime_root

        self._machobj = Machines(machine=machine_name)

        self._no_setup = no_setup
        self._no_build = no_build or no_setup or namelists_only
        self._no_run   = no_run or self._no_build

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
        self._test_root = self._machobj.get_value("CESMSCRATCHROOT") if test_root is None else test_root
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
                else self._machobj.get_value("CCSM_BASELINE")

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
        # Since the name-list phase can fail without aborting later phases, we
        # need some extra state to remember tests that had namelist problems.
        # name -> (phase, status, has_namelist_problem)
        self._tests = {}
        for test_name in test_names:
            self._tests[test_name] = (INITIAL_PHASE, TEST_PASS_STATUS, False)

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
        if not self._baseline_cmp_name and not self._baseline_gen_name:
            self._phases.remove(NAMELIST_PHASE)

        if use_existing:
            for test in self._tests:
                ts = TestStatus(self._get_test_dir(test))
                for phase, status in ts:
                    self._update_test_status(test, phase, TEST_PENDING_STATUS)
                    self._update_test_status(test, phase, status)
        else:
            # None of the test directories should already exist.
            for test in self._tests:
                expect(not os.path.exists(self._get_test_dir(test)),
                       "Cannot create new case in directory '%s', it already exists."
                       " Pick a different test-id" % self._get_test_dir(test))

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
        return status not in CONTINUE and status != TEST_PENDING_STATUS

    ###########################################################################
    def _work_remains(self, test):
    ###########################################################################
        test_phase, test_status, _ = self._get_test_data(test)
        return (test_status in CONTINUE or test_status == TEST_PENDING_STATUS) and\
            test_phase != self._phases[-1]

    ###########################################################################
    def _get_test_status(self, test, phase=None):
    ###########################################################################
        curr_phase, curr_status, nl_fail = self._get_test_data(test)
        if phase == NAMELIST_PHASE and nl_fail:
            return NAMELIST_FAIL_STATUS
        elif phase is None or phase == curr_phase:
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
        old_phase, old_status, old_nl_fail = self._get_test_data(test)

        if old_phase == phase:
            expect(old_status == TEST_PENDING_STATUS,
                   "Only valid to transition from PENDING to something else, found '%s' for phase '%s'" %
                   (old_status, phase))
            expect(status != TEST_PENDING_STATUS,
                   "Cannot transition from PEND -> PEND")
        else:
            expect(old_status in CONTINUE,
                   "Why did we move on to next phase when prior phase did not pass?")
            expect(status == TEST_PENDING_STATUS,
                   "New phase should be set to pending status")
            expect(self._phases.index(old_phase) == phase_idx - 1,
                   "Skipped phase?")
        # Must be atomic
        self._tests[test] = (phase, status, old_nl_fail)

    ###########################################################################
    def _test_has_nl_problem(self, test):
    ###########################################################################
        curr_phase, curr_status, _ = self._get_test_data(test)
        expect(curr_phase == NAMELIST_PHASE, "Setting namelist status outside of namelist phase?")
        # Must be atomic
        self._tests[test] = (curr_phase, curr_status, True)

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
                                     "%s PASSED for test '%s'.\nCommand: %s\nOutput: %s\n\nErrput: %s" %
                                     (phase, test, cmd, output, errput))
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
        elif test in self._test_data and "wallclock" in self._test_data[test]:
            create_newcase_cmd += " --walltime %s" % self._test_data[test]['wallclock']

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
        drv_config_file = files.get_value("CONFIG_DRV_FILE")
        drv_comp = Component(drv_config_file)
        envtest.add_elements_by_group(files, {}, "env_test.xml")
        envtest.add_elements_by_group(drv_comp, {}, "env_test.xml")
        envtest.set_value("TESTCASE", test_case)
        envtest.set_value("TEST_TESTID", self._test_id)
        envtest.set_value("CASEBASEID", test)

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
        lockedfiles = os.path.join(test_dir, "LockedFiles")
        if not os.path.exists(lockedfiles):
            os.mkdir(lockedfiles)
        shutil.copy(os.path.join(test_dir,"env_run.xml"),
                    os.path.join(lockedfiles, "env_run.orig.xml"))

        with Case(test_dir, read_only=False) as case:
            case.set_value("SHAREDLIBROOT",
                           os.path.join(self._test_root,
                                        "sharedlibroot.%s"%self._test_id))
            envtest.set_initial_values(case)
            if self._save_timing:
                case.set_value("SAVE_TIMING", True)

        return True

    ###########################################################################
    def _setup_phase(self, test):
    ###########################################################################
        test_dir  = self._get_test_dir(test)

        return self._shell_cmd_for_phase(test, "./case.setup", SETUP_PHASE, from_dir=test_dir)

    ###########################################################################
    def _nlcomp_phase(self, test):
    ###########################################################################
        test_dir       = self._get_test_dir(test)
        casedoc_dir    = os.path.join(test_dir, "CaseDocs")
        compare_nl     = os.path.join(CIME.utils.get_scripts_root(), "Tools", "compare_namelists")
        simple_compare = os.path.join(CIME.utils.get_scripts_root(), "Tools", "simple_compare")

        if self._baseline_cmp_name:
            has_fails         = False
            baseline_dir      = os.path.join(self._baseline_root, self._baseline_cmp_name, test)
            baseline_casedocs = os.path.join(baseline_dir, "CaseDocs")

            # Start off by comparing everything in CaseDocs except a few arbitrary files (ugh!)
            # TODO: Namelist files should have consistent suffix
            all_items_to_compare = [item for item in glob.glob("%s/*" % casedoc_dir)\
                                    if "README" not in os.path.basename(item)\
                                    and not item.endswith("doc")\
                                    and not item.endswith("prescribed")\
                                    and not os.path.basename(item).startswith(".")] + \
                                    glob.glob("%s/*user_nl*" % test_dir)
            for item in all_items_to_compare:
                baseline_counterpart = os.path.join(baseline_casedocs \
                                                    if os.path.dirname(item).endswith("CaseDocs") \
                                                    else baseline_dir,os.path.basename(item))
                if not os.path.exists(baseline_counterpart):
                    self._log_output(test, "Missing baseline namelist '%s'" % baseline_counterpart)
                    has_fails = True
                else:
                    if CIME.compare_namelists.is_namelist_file(item):
                        rc, output, _  = run_cmd("%s %s %s -c %s 2>&1" %
                                                 (compare_nl, baseline_counterpart, item, test))
                    else:
                        rc, output, _  = run_cmd("%s %s %s -c %s 2>&1" %
                                                 (simple_compare, baseline_counterpart, item, test))

                    if rc != 0:
                        has_fails = True
                        self._log_output(test, output)

            if has_fails:
                self._test_has_nl_problem(test)

        if self._baseline_gen_name:
            baseline_dir      = os.path.join(self._baseline_root, self._baseline_gen_name, test)
            baseline_casedocs = os.path.join(baseline_dir, "CaseDocs")
            if not os.path.isdir(baseline_dir):
                os.makedirs(baseline_dir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH)

            if os.path.isdir(baseline_casedocs):
                shutil.rmtree(baseline_casedocs)

            shutil.copytree(casedoc_dir, baseline_casedocs)
            os.chmod(baseline_casedocs, stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH)
            for item in glob.glob("%s/*" % baseline_casedocs):
                os.chmod(item, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP)

            for item in glob.glob(os.path.join(test_dir, "user_nl*")):
                preexisting_baseline = os.path.join(baseline_dir, os.path.basename(item))
                if (os.path.exists(preexisting_baseline)):
                    os.remove(preexisting_baseline)
                shutil.copy2(item, baseline_dir)
                os.chmod(preexisting_baseline, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP)

        # Always mark as passed unless we hit exception
        return True

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
    def _get_procs_needed(self, test, phase, threads_in_flight=None):
    ###########################################################################
        if phase == RUN_PHASE and self._no_batch:
            test_dir = self._get_test_dir(test)
            out = run_cmd_no_fail("./xmlquery TOTALPES -value", from_dir=test_dir)
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
        status  = (TEST_PENDING_STATUS if test_phase == RUN_PHASE and not \
                   self._no_batch else TEST_PASS_STATUS) if success else TEST_FAIL_STATUS

        if status != TEST_PENDING_STATUS:
            self._update_test_status(test, test_phase, status)

        status_str = "Finished %s for test %s in %f seconds (%s)" %\
                     (test_phase, test, elapsed_time, status)
        if not success:
            status_str += "    Case dir: %s" % self._get_test_dir(test)
        logger.info(status_str)

        if test_phase in [CREATE_NEWCASE_PHASE, XML_PHASE, NAMELIST_PHASE]:
            # These are the phases for which TestScheduler is reponsible for
            # updating the TestStatus file
            nl_problem = self._get_test_data(test)[2]
            status = TEST_FAIL_STATUS if nl_problem and test_phase == NAMELIST_PHASE else status
            self._update_test_status_file(test, test_phase, status)

        # On batch systems, we want to immediately submit to the queue, because
        # it's very cheap to submit and will get us a better spot in line
        if (success and not self._no_run and not self._no_batch and test_phase == MODEL_BUILD_PHASE):
            logger.info("Starting %s for test %s with %d procs" % (RUN_PHASE, test, 1))
            self._update_test_status(test, RUN_PHASE, TEST_PENDING_STATUS)
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
                        test_phase, test_status, _ = self._get_test_data(test)
                        expect(test_status != TEST_PENDING_STATUS, test)
                        next_phase = self._phases[self._phases.index(test_phase) + 1]
                        procs_needed = self._get_procs_needed(test, next_phase, threads_in_flight)

                        if procs_needed <= self._procs_avail:
                            self._procs_avail -= procs_needed

                            # Necessary to print this way when multiple threads printing
                            logger.info("Starting %s for test %s with %d procs" %
                                        (next_phase, test, procs_needed))

                            self._update_test_status(test, next_phase, TEST_PENDING_STATUS)
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
                                self._update_test_status(test, next_phase, TEST_PENDING_STATUS)
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
    def run_tests(self):
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

        # TODO - documentation

        self._producer()

        expect(threading.active_count() == 1, "Leftover threads?")

        # Setup cs files
        self._setup_cs_files()

        # Return True if all tests passed from our point of view
        logger.info( "At test-scheduler close, state is:")
        rv = True
        for test in self._tests:
            phase, status, nl_fail = self._get_test_data(test)

            if status == TEST_PASS_STATUS and phase == RUN_PHASE:
                # Be cautious about telling the user that the test passed. This
                # status should match what they would see on the dashboard. Our
                # self._test_states does not include comparison fail information,
                # so we need to parse test status.
                ts = TestStatus(self._get_test_dir(test))
                status = ts.get_overall_test_status()

            if status not in [TEST_PASS_STATUS, TEST_PENDING_STATUS]:
                logger.info( "%s %s (phase %s)" % (status, test, phase))
                rv = False

            elif nl_fail:
                logger.info( "%s %s (but otherwise OK)" % (NAMELIST_FAIL_STATUS, test))
                rv = False

            else:
                logger.info("%s %s %s" % (status, test, phase))

            logger.info( "    Case dir: %s" % self._get_test_dir(test))

        logger.info( "test-scheduler took %s seconds"% (time.time() - start_time))

        return rv
