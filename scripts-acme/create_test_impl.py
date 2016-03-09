"""
Implementation of create_test functionality from CIME
"""

import sys, os, shutil, traceback, stat, glob, threading, time, thread

import acme_util, compare_namelists, wait_for_tests

from acme_util import expect, warning, verbose_print, run_cmd
from wait_for_tests import TEST_PASS_STATUS, TEST_FAIL_STATUS, TEST_PENDING_STATUS, TEST_STATUS_FILENAME, NAMELIST_FAIL_STATUS, RUN_PHASE, NAMELIST_PHASE

INITIAL_PHASE = "INIT"
CREATE_NEWCASE_PHASE = "CREATE_NEWCASE"
XML_PHASE   = "XML"
SETUP_PHASE = "SETUP"
BUILD_PHASE = "BUILD"
TEST_STATUS_PHASE = "TEST_STATUS"
PHASES = [INITIAL_PHASE, CREATE_NEWCASE_PHASE, XML_PHASE, SETUP_PHASE, NAMELIST_PHASE, BUILD_PHASE, RUN_PHASE] # Order matters
CONTINUE = [TEST_PASS_STATUS, NAMELIST_FAIL_STATUS]

###############################################################################
class CreateTest(object):
###############################################################################

    ###########################################################################
    def __init__(self, test_names,
                 no_run=False, no_build=False, no_batch=None,
                 test_root=None, test_id=None,
                 compiler=None,
                 baseline_root=None, baseline_name=None,
                 clean=False,
                 compare=False, generate=False, namelists_only=False,
                 project=None, parallel_jobs=None):
    ###########################################################################
        self._cime_root      = acme_util.get_cime_root()
        self._test_names     = test_names
        self._no_build       = no_build      if not namelists_only else True
        self._no_run         = no_run        if not self._no_build else True
        self._no_batch       = no_batch      if no_batch is not None else not acme_util.does_machine_have_batch()
        self._test_root      = test_root     if test_root is not None else acme_util.get_machine_info("CESMSCRATCHROOT")
        self._test_id        = test_id       if test_id is not None else acme_util.get_utc_timestamp()
        self._project        = project       if project is not None else acme_util.get_machine_project()
        self._baseline_root  = baseline_root if baseline_root is not None else acme_util.get_machine_info("CCSM_BASELINE", project=self._project)
        self._baseline_name  = None
        self._compiler       = compiler      if compiler is not None else acme_util.get_machine_info("COMPILERS")[0]
        self._clean          = clean
        self._compare        = compare
        self._generate       = generate
        self._namelists_only = namelists_only
        self._parallel_jobs  = parallel_jobs if parallel_jobs is not None else min(len(self._test_names), int(acme_util.get_machine_info("MAX_TASKS_PER_NODE")))

        # Oversubscribe by 1/4
        pes = int(acme_util.get_machine_info("MAX_TASKS_PER_NODE"))

        # This is the only data that multiple threads will simultaneously access
        # Each test has it's own index and setting/retrieving items from a list
        # is atomic, so this should be fine to use without mutex
        self._test_states    = [ (INITIAL_PHASE, TEST_PASS_STATUS) ] * len(test_names)

        self._proc_pool = int(pes * 1.25)

        # Since the name-list phase can fail without aborting later phases, we
        # need some extra state to remember tests that had namelist problems
        self._tests_with_nl_problems = [None] * len(test_names)

        # Setup phases
        self._phases = list(PHASES)
        if (no_build):
            self._phases.remove(BUILD_PHASE)
        if (no_run):
            self._phases.remove(RUN_PHASE)

        if (not self._compare and not self._generate):
            self._phases.remove(NAMELIST_PHASE)
        else:
            if (baseline_name is None):
                branch_name = acme_util.get_current_branch(repo=self._cime_root)
                expect(branch_name is not None, "Could not determine baseline name from branch, please use -b option")
                self._baseline_name = os.path.join(self._compiler, branch_name)
            else:
                self._baseline_name  = baseline_name
                if (not self._baseline_name.startswith("%s/" % self._compiler)):
                    self._baseline_name = os.path.join(self._compiler, self._baseline_name)

        # Validate any assumptions that were not caught by the arg parser

        # None of the test directories should already exist.
        for test in self._test_names:
            expect(not os.path.exists(self._get_test_dir(test)),
                   "Cannot create new case in directory '%s', it already exists. Pick a different test-id" % self._get_test_dir(test))

        # By the end of this constructor, this program should never hard abort,
        # instead, errors will be placed in the TestStatus files for the various
        # tests cases


    ###########################################################################
    def _log_output(self, test_name, output):
    ###########################################################################
        test_dir = self._get_test_dir(test_name)
        if (not os.path.isdir(test_dir)):
            # Note: making this directory could cause create_newcase to fail
            # if this is run before.
            os.makedirs(test_dir)

        with open(os.path.join(test_dir, "TestStatus.log"), 'a') as fd:
            fd.write(output)

    ###########################################################################
    def _get_case_id(self, test_name):
    ###########################################################################
        baseline_action_code = ".C" if self._compare else (".G" if self._generate else "")
        return "%s%s.%s" % (test_name, baseline_action_code, self._test_id)

    ###########################################################################
    def _get_test_dir(self, test_name):
    ###########################################################################
        return os.path.join(self._test_root, self._get_case_id(test_name))

    ###########################################################################
    def _get_test_data(self, test_name):
    ###########################################################################
        state_idx = self._test_names.index(test_name)
        return self._test_states[state_idx]

    ###########################################################################
    def _is_broken(self, test_name):
    ###########################################################################
        status = self._get_test_status(test_name)
        return status not in CONTINUE and status != TEST_PENDING_STATUS

    ###########################################################################
    def _work_remains(self, test_name):
    ###########################################################################
        test_phase, test_status = self._get_test_data(test_name)
        return (test_status in CONTINUE or test_status == TEST_PENDING_STATUS) and test_phase != self._phases[-1]

    ###########################################################################
    def _get_test_status(self, test_name, phase=None):
    ###########################################################################
        curr_phase = self._get_test_phase(test_name)
        if (phase == NAMELIST_PHASE and test_name in self._tests_with_nl_problems):
            return NAMELIST_FAIL_STATUS
        elif (phase is None or phase == curr_phase):
            return self._get_test_data(test_name)[1]
        else:
            expect(phase is None or self._phases.index(phase) < self._phases.index(curr_phase),
                   "Tried to see the future")
            # Assume all older phases PASSed
            return TEST_PASS_STATUS

    ###########################################################################
    def _get_test_phase(self, test_name):
    ###########################################################################
        return self._get_test_data(test_name)[0]

    ###########################################################################
    def _update_test_status(self, test_name, phase, status):
    ###########################################################################
        state_idx = self._test_names.index(test_name)
        phase_idx = self._phases.index(phase)
        old_phase, old_status = self._test_states[state_idx]

        if (old_phase == phase):
            expect(old_status == TEST_PENDING_STATUS,
                   "Only valid to transition from PENDING to something else, found '%s'" % old_status)
            expect(status != TEST_PENDING_STATUS,
                   "Cannot transition from PEND -> PEND")
        else:
            expect(old_status in CONTINUE,
                   "Why did we move on to next phase when prior phase did not pass?")
            expect(status == TEST_PENDING_STATUS,
                   "New phase should be set to pending status")
            expect(self._phases.index(old_phase) == phase_idx - 1,
                   "Skipped phase?")

        self._test_states[state_idx] = (phase, status)

    ###########################################################################
    def _run_phase_command(self, test_name, cmd, phase, from_dir=None):
    ###########################################################################
        while (True):
            rc, output, errput = run_cmd(cmd, ok_to_fail=True, from_dir=from_dir)
            if (rc != 0):
                self._log_output(test_name,
                                 "%s FAILED for test '%s'.\nCommand: %s\nOutput: %s\n\nErrput: %s" %
                                 (phase, test_name, cmd, output, errput))
                # Temporary hack to get around odd file descriptor use by
                # buildnml scripts.
                if ("bad interpreter" in errput):
                    time.sleep(1)
                    continue
                else:
                    break
            else:
                self._log_output(test_name,
                                 "%s PASSED for test '%s'.\nCommand: %s\nOutput: %s\n\nErrput: %s" %
                                 (phase, test_name, cmd, output, errput))
                break

        return rc == 0

    ###########################################################################
    def _create_newcase_phase(self, test_name):
    ###########################################################################
        test_dir = self._get_test_dir(test_name)

        test_case, case_opts, grid, compset, machine, compiler, test_mods = acme_util.parse_test_name(test_name)
        if (compiler != self._compiler):
            raise StandardError("Test '%s' has compiler that does not match instance compliler '%s'" % (test_name, self._compiler))
        if (self._parallel_jobs == 1):
            scratch_dir = acme_util.get_machine_info("CESMSCRATCHROOT", machine=machine, project=self._project)
            sharedlibroot = os.path.join(scratch_dir, "sharedlibroot.%s" % self._test_id)
        else:
            # Parallelizing builds introduces potential sync problems with sharedlibroot
            # Just let every case build it's own
            sharedlibroot = os.path.join(test_dir, "sharedlibroot.%s" % self._test_id)

        create_newcase_cmd = "%s -silent -case %s -res %s -mach %s -compiler %s -compset %s -testname %s -project %s -nosavetiming -sharedlibroot %s" % \
                              (os.path.join(self._cime_root, "scripts", "create_newcase"),
                               test_dir, grid, machine, compiler, compset, test_case, self._project,
                               sharedlibroot)
        if (case_opts is not None):
            create_newcase_cmd += " -confopts _%s" % ("_".join(case_opts))
        if (test_mods is not None):
            test_mod_file = os.path.join(self._cime_root, "scripts", "Testing", "Testlistxml", "testmods_dirs", test_mods)
            if (not os.path.exists(test_mod_file)):
                self._log_output(test_name, "Missing testmod file '%s'" % test_mod_file)
                return False
            create_newcase_cmd += " -user_mods_dir %s" % test_mod_file

        return self._run_phase_command(test_name, create_newcase_cmd, CREATE_NEWCASE_PHASE)

    ###########################################################################
    def _xml_phase(self, test_name):
    ###########################################################################
        test_case, _, _, _, machine, _, _ = acme_util.parse_test_name(test_name)

        xml_file = os.path.join(self._get_test_dir(test_name), "env_test.xml")
        xml_bridge_cmd = os.path.join(acme_util.get_acme_scripts_root(), "xml_bridge")

        mach_dir = os.path.join(self._cime_root, "machines-acme")
        xml_bridge_cmd += " %s %s %s" % (mach_dir, machine, xml_file)

        xml_bridge_cmd += " TESTCASE,%s" % test_case
        xml_bridge_cmd += " TEST_TESTID,%s" % self._test_id

        test_argv = "-testname %s -testroot %s" % (test_name, self._test_root)
        if (self._generate):
            test_argv += " -generate %s" % self._baseline_name
        if (self._compare):
            test_argv += " -compare %s" % self._baseline_name
        xml_bridge_cmd += " 'TEST_ARGV,%s'" % test_argv

        xml_bridge_cmd += " CASEBASEID,%s" % test_name

        if (self._generate):
            xml_bridge_cmd += " BASELINE_NAME_GEN,%s" % self._baseline_name
            xml_bridge_cmd += " BASEGEN_CASE,%s" % os.path.join(self._baseline_name, test_name)
        if (self._compare):
            xml_bridge_cmd += " BASELINE_NAME_CMP,%s" % self._baseline_name
            xml_bridge_cmd += " BASECMP_CASE,%s" % os.path.join(self._baseline_name, test_name)

        xml_bridge_cmd += " CLEANUP,%s" % ("TRUE" if self._clean else "FALSE")

        if (self._generate or self._compare):
            xml_bridge_cmd += " BASELINE_ROOT,%s" % self._baseline_root

        xml_bridge_cmd += " GENERATE_BASELINE,%s" % ("TRUE" if self._generate else "FALSE")
        xml_bridge_cmd += " COMPARE_BASELINE,%s" % ("TRUE" if self._compare else "FALSE")

        return self._run_phase_command(test_name, xml_bridge_cmd, XML_PHASE)

    ###########################################################################
    def _setup_phase(self, test_name):
    ###########################################################################
        test_case = acme_util.parse_test_name(test_name)[0]
        test_dir  = self._get_test_dir(test_name)
        test_case_definition_dir = os.path.join(self._cime_root, "scripts", "Testing", "Testcases")
        test_build = os.path.join(test_dir, "%s.test_build" % self._get_case_id(test_name))

        if (os.path.exists(os.path.join(test_case_definition_dir, "%s_build.csh" % test_case))):
            shutil.copy(os.path.join(test_case_definition_dir, "%s_build.csh" % test_case), test_build)
        else:
            shutil.copy(os.path.join(test_case_definition_dir, "tests_build.csh"), test_build)

        return self._run_phase_command(test_name, "./cesm_setup", SETUP_PHASE, from_dir=test_dir)

    ###########################################################################
    def _nlcomp_phase(self, test_name):
    ###########################################################################
        test_dir          = self._get_test_dir(test_name)
        casedoc_dir       = os.path.join(test_dir, "CaseDocs")
        baseline_dir      = os.path.join(self._baseline_root, self._baseline_name, test_name)
        baseline_casedocs = os.path.join(baseline_dir, "CaseDocs")
        compare_nl        = os.path.join(test_dir, "Tools/compare_namelists")
        simple_compare    = os.path.join(test_dir, "Tools/simple_compare")

        if (self._compare):
            has_fails = False

            # Start off by comparing everything in CaseDocs except a few arbitrary files (ugh!)
            # TODO: Namelist files should have consistent suffix
            all_items_to_compare = \
                [ item for item in glob.glob("%s/*" % casedoc_dir) if "README" not in os.path.basename(item) and not item.endswith("doc") and not item.endswith("prescribed") and not os.path.basename(item).startswith(".")] + \
                glob.glob("%s/*user_nl*" % test_dir)
            for item in all_items_to_compare:
                baseline_counterpart = os.path.join(baseline_casedocs if os.path.dirname(item).endswith("CaseDocs") else baseline_dir,
                                                    os.path.basename(item))
                if (not os.path.exists(baseline_counterpart)):
                    self._log_output(test_name, "Missing baseline namelist '%s'" % baseline_counterpart)
                    has_fails = True
                else:
                    if (compare_namelists.is_namelist_file(item)):
                        rc, output, _  = run_cmd("%s %s %s -c %s 2>&1" % (compare_nl, baseline_counterpart, item, test_name), ok_to_fail=True)
                    else:
                        rc, output, _  = run_cmd("%s %s %s -c %s 2>&1" % (simple_compare, baseline_counterpart, item, test_name), ok_to_fail=True)

                    if (rc != 0):
                        has_fails = True
                        self._log_output(test_name, output)

            if (has_fails):
                idx = self._test_names.index(test_name)
                self._tests_with_nl_problems[idx] = test_name

        elif (self._generate):
            if (not os.path.isdir(baseline_dir)):
                os.makedirs(baseline_dir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH)

            if (os.path.isdir(baseline_casedocs)):
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
    def _build_phase(self, test_name):
    ###########################################################################
        case_id = self._get_case_id(test_name)
        test_dir = self._get_test_dir(test_name)
        return self._run_phase_command(test_name, "./%s.test_build" % case_id, BUILD_PHASE, from_dir=test_dir)

    ###########################################################################
    def _run_phase(self, test_name):
    ###########################################################################
        case_id = self._get_case_id(test_name)
        test_dir = self._get_test_dir(test_name)
        if (self._no_batch):
            return self._run_phase_command(test_name, "./%s.test" % case_id, RUN_PHASE, from_dir=test_dir)
        else:
            return self._run_phase_command(test_name, "./%s.submit" % case_id, RUN_PHASE, from_dir=test_dir)

    ###########################################################################
    def _update_test_status_file(self, test_name):
    ###########################################################################
        # TODO: The run scripts heavily use the TestStatus file. So we write out
        # the phases we have taken care of and then let the run scrips go from there
        # Eventually, it would be nice to have TestStatus management encapsulated
        # into a single place.

        str_to_write = ""
        made_it_to_phase = self._get_test_phase(test_name)
        made_it_to_phase_idx = self._phases.index(made_it_to_phase)
        for phase in self._phases[0:made_it_to_phase_idx+1]:
            str_to_write += "%s %s %s\n" % (self._get_test_status(test_name, phase), test_name, phase)

        if (not self._no_run and not self._is_broken(test_name) and made_it_to_phase == BUILD_PHASE):
            # Ensure PEND state always gets added to TestStatus file if we are
            # about to run test
            str_to_write += "%s %s %s\n" % (TEST_PENDING_STATUS, test_name, RUN_PHASE)

        test_status_file = os.path.join(self._get_test_dir(test_name), TEST_STATUS_FILENAME)
        with open(test_status_file, "w") as fd:
            fd.write(str_to_write)

    ###########################################################################
    def _run_catch_exceptions(self, test_name, phase, run):
    ###########################################################################
        try:
            return run(test_name)
        except Exception as e:
            exc_tb = sys.exc_info()[2]
            errput = "Test '%s' failed in phase '%s' with exception '%s'" % (test_name, phase, str(e))
            self._log_output(test_name, errput)
            warning("Caught exception: %s" % str(e))
            traceback.print_tb(exc_tb)
            return False

    ###########################################################################
    def _get_procs_needed(self, test_name, phase):
    ###########################################################################
        if (phase == RUN_PHASE and self._no_batch):
            test_dir = self._get_test_dir(test_name)
            out = run_cmd("./xmlquery TOTALPES -value", from_dir=test_dir)
            return int(out)
        else:
            return 1

    ###########################################################################
    def _handle_test_status_file(self, test_name, test_phase, success):
    ###########################################################################
        #
        # This complexity is due to sharing of TestStatus responsibilities
        #

        try:
            if (test_phase != RUN_PHASE and
                (not success or test_phase == BUILD_PHASE or test_phase == self._phases[-1])):
                self._update_test_status_file(test_name)

            # If we failed VERY early on in the run phase, it's possible that
            # the CIME scripts never got a chance to set the state.
            elif (test_phase == RUN_PHASE and not success):
                test_status_file = os.path.join(self._get_test_dir(test_name), TEST_STATUS_FILENAME)

                statuses = wait_for_tests.parse_test_status_file(test_status_file)[0]
                if ( RUN_PHASE not in statuses or
                     statuses[RUN_PHASE] in [TEST_PASS_STATUS, TEST_PENDING_STATUS] ):
                    self._update_test_status_file(test_name)

        except Exception as e:
            # TODO: What to do here? This failure is very severe because the
            # only way for test results to be communicated is by the TestStatus
            # file.
            warning("VERY BAD! Could not handle TestStatus file '%s': '%s'" %
                    (os.path.join(self._get_test_dir(test_name), TEST_STATUS_FILENAME), str(e)))
            thread.interrupt_main()

    ###########################################################################
    def _wait_for_something_to_finish(self, threads_in_flight):
    ###########################################################################
        expect(len(threads_in_flight) <= self._parallel_jobs, "Oversubscribed?")
        finished_tests = []
        while (not finished_tests):
            for test_name, thread_info in threads_in_flight.iteritems():
                if (not thread_info[0].is_alive()):
                    finished_tests.append( (test_name, thread_info[1]) )

            if (not finished_tests):
                time.sleep(0.2)

        for finished_test, procs_needed in finished_tests:
            self._proc_pool += procs_needed
            del threads_in_flight[finished_test]

    ###########################################################################
    def _consumer(self, test_name, test_phase, phase_method):
    ###########################################################################
        before_time = time.time()
        success = self._run_catch_exceptions(test_name, test_phase, phase_method)
        elapsed_time = time.time() - before_time
        status  = (TEST_PENDING_STATUS if test_phase == RUN_PHASE and not self._no_batch else TEST_PASS_STATUS) if success else TEST_FAIL_STATUS

        if (status != TEST_PENDING_STATUS):
            self._update_test_status(test_name, test_phase, status)
        self._handle_test_status_file(test_name, test_phase, success)

        status_str = "Finished %s for test %s in %f seconds (%s)\n" % (test_phase, test_name, elapsed_time, status)
        if (not success):
            status_str += "    Case dir: %s\n" % self._get_test_dir(test_name)
        sys.stdout.write(status_str)

    ###########################################################################
    def _producer(self):
    ###########################################################################
        threads_in_flight = {} # test-name -> (thread, procs)
        while (True):
            work_to_do = False
            num_threads_launched_this_iteration = 0
            for test_name in self._test_names:

                # If we have no workers available, immediately wait
                if (len(threads_in_flight) == self._parallel_jobs):
                    self._wait_for_something_to_finish(threads_in_flight)

                if (self._work_remains(test_name)):
                    work_to_do = True
                    if (test_name not in threads_in_flight):
                        test_phase, test_status = self._get_test_data(test_name)
                        expect(test_status != TEST_PENDING_STATUS, test_name)
                        next_phase = self._phases[self._phases.index(test_phase) + 1]
                        procs_needed = self._get_procs_needed(test_name, next_phase)

                        if (procs_needed <= self._proc_pool):
                            self._proc_pool -= procs_needed

                            # Necessary to print this way when multiple threads printing
                            sys.stdout.write("Starting %s for test %s with %d procs\n" % (next_phase, test_name, procs_needed))

                            self._update_test_status(test_name, next_phase, TEST_PENDING_STATUS)
                            t = threading.Thread(target=self._consumer,
                                                 args=(test_name, next_phase, getattr(self, "_%s_phase" % next_phase.lower()) ))
                            threads_in_flight[test_name] = (t, procs_needed)
                            t.start()
                            num_threads_launched_this_iteration += 1

            if (not work_to_do):
                break

            if (num_threads_launched_this_iteration == 0):
                # No free resources, wait for something in flight to finish
                self._wait_for_something_to_finish(threads_in_flight)

        for thread_info in threads_in_flight.values():
            thread_info[0].join()

    ###########################################################################
    def _setup_cs_files(self):
    ###########################################################################
        try:
            acme_scripts_root = acme_util.get_acme_scripts_root()
            template_file = os.path.join(acme_scripts_root, "cs.status.template")
            template = open(template_file, "r").read()
            template = template.replace("<PATH>", acme_scripts_root).replace("<TESTID>", self._test_id)

            cs_status_file = os.path.join(self._test_root, "cs.status.%s" % self._test_id)
            with open(cs_status_file, "w") as fd:
                fd.write(template)
            os.chmod(cs_status_file, os.stat(cs_status_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)

            template_file = os.path.join(acme_scripts_root, "cs.submit.template")
            template = open(template_file, "r").read()
            build_cmd = "./*.test_build" if self._no_build else ":"
            run_cmd = "./*.test" if self._no_batch else "./*.submit"
            template = template.replace("<BUILD_CMD>", build_cmd).replace("<RUN_CMD>", run_cmd).replace("<TESTID>", self._test_id)

            if (self._no_build or self._no_run):
                cs_submit_file = os.path.join(self._test_root, "cs.submit.%s" % self._test_id)
                with open(cs_submit_file, "w") as fd:
                    fd.write(template)
                os.chmod(cs_submit_file, os.stat(cs_submit_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)

        except Exception as e:
            warning("FAILED to set up cs files: %s" % str(e))

    ###########################################################################
    def create_test(self):
    ###########################################################################
        """
        Main API for this class.

        Return True if all tests passed.
        """
        start_time = time.time()

        # Tell user what will be run
        print "RUNNING TESTS:"
        for test_name in self._test_names:
            print " ", test_name

        # TODO - documentation

        self._producer()

        expect(threading.active_count() == 1, "Leftover threads?")

        # Setup cs files
        self._setup_cs_files()

        # Return True if all tests passed
        print "At create_test close, state is:"
        rv = True
        for idx, test_name in enumerate(self._test_names):
            phase, status = self._test_states[idx]

            if (status == TEST_PASS_STATUS and phase == RUN_PHASE):
                # Be cautious about telling the user that the test passed. This
                # status should match what they would see on the dashboard. Our
                # self._test_states does not include comparison fail information,
                # so we need to parse test status.
                test_status_file = os.path.join(self._get_test_dir(test_name), TEST_STATUS_FILENAME)
                status = wait_for_tests.interpret_status_file(test_status_file)[1]

            if (status not in [TEST_PASS_STATUS, TEST_PENDING_STATUS]):
                print "%s %s (phase %s)" % (status, test_name, phase)
                rv = False

            elif (test_name in self._tests_with_nl_problems):
                print "%s %s (but otherwise OK)" % (NAMELIST_FAIL_STATUS, test_name)
                rv = False

            else:
                print status, test_name, phase

            print "    Case dir: %s" % self._get_test_dir(test_name)

        print "create_test took", time.time() - start_time, "seconds"

        return rv
