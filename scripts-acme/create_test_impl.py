"""
Implementation of create_test functionality from CIME
"""

import sys, os, shutil, traceback, stat, glob

import acme_util, compare_namelists, wait_for_tests

from acme_util import expect, warning, verbose_print, run_cmd
from wait_for_tests import TEST_PASSED_STATUS, TEST_FAIL_STATUS, TEST_PENDING_STATUS, TEST_STATUS_FILENAME, NAMELIST_FAIL_STATUS, RUN_PHASE, NAMELIST_PHASE

INITIAL_PHASE = "INITIAL_PHASE"
CREATE_NEWCASE_PHASE = "CREATE_NEWCASE"
XML_PHASE   = "XML_PHASE"
SETUP_PHASE = "SETUP"
BUILD_PHASE = "BUILD"
PHASES = [INITIAL_PHASE, CREATE_NEWCASE_PHASE, XML_PHASE, SETUP_PHASE, NAMELIST_PHASE, BUILD_PHASE, RUN_PHASE] # Order matters

###############################################################################
class CreateTest(object):
###############################################################################

    ###########################################################################
    def __init__(self, test_names,
                 no_run, no_build, no_batch,
                 test_root, test_id,
                 baseline_root, baseline_name,
                 clean,
                 compare, generate, namelists_only,
                 project):
    ###########################################################################
        self._test_names     = test_names
        self._no_run         = no_run
        self._no_build       = no_build
        self._no_batch       = no_batch
        self._test_root      = test_root
        self._test_id        = test_id
        self._baseline_root  = baseline_root
        self._baseline_name  = baseline_name
        self._clean          = clean
        self._compare        = compare
        self._generate       = generate
        self._namelists_only = namelists_only
        self._project        = project

        self._cime_root      = acme_util.get_cime_root()

        self._test_states    = [ (INITIAL_PHASE, TEST_PASSED_STATUS) ] * len(test_names)

        # Since the name-list phase can fail without aborting later phases, we
        # need some extra state to remember tests that had namelist problems
        self._tests_with_nl_problems = []

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
    def _get_test_data(self, test_name, idx):
    ###########################################################################
        state_idx = self._test_names.index(test_name)
        return self._test_states[state_idx][idx]

    ###########################################################################
    def _is_broken(self, test_name):
    ###########################################################################
        return self._get_test_status(test_name) not in [TEST_PASSED_STATUS, NAMELIST_FAIL_STATUS, TEST_PENDING_STATUS]

    ###########################################################################
    def _get_test_status(self, test_name, phase=None):
    ###########################################################################
        if (phase == NAMELIST_PHASE and test_name in self._tests_with_nl_problems):
            return NAMELIST_FAIL_STATUS
        elif (phase is None or phase == self._get_test_phase(test_name)):
            return self._get_test_data(test_name, 1)
        else:
            expect(phase is None or PHASES.index(phase) < PHASES.index(self._get_test_phase(test_name)),
                   "Tried to see the future")
            # Assume all older phases PASSed
            return TEST_PASSED_STATUS

    ###########################################################################
    def _get_test_phase(self, test_name):
    ###########################################################################
        return self._get_test_data(test_name, 0)

    ###########################################################################
    def _update_test_status(self, test_name, phase, status):
    ###########################################################################
        state_idx = self._test_names.index(test_name)
        phase_idx = PHASES.index(phase)

        expect(not self._is_broken(test_name),
               "Why did we move on to next phase when prior phase did not pass?")
        expect(PHASES.index(self._get_test_phase(test_name)) == phase_idx - 1,
               "Skipped phase?")

        self._test_states[state_idx] = (phase, status)

    ###########################################################################
    def _run_phase_command(self, test_name, cmd, phase, from_dir=None, pass_status=TEST_PASSED_STATUS):
    ###########################################################################
        rc, output, errput = run_cmd(cmd, ok_to_fail=True, from_dir=from_dir)
        if (rc != 0):
            self._log_output(test_name,
                             "%s FAILED for test '%s'.\nCommand: %s\nOutput: %s\n\nErrput: %s" %
                             (phase, test_name, cmd, output, errput))
        else:
            self._log_output(test_name,
                             "%s PASSED for test '%s'.\nCommand: %s\nOutput: %s\n\nErrput: %s" %
                             (phase, test_name, cmd, output, errput))

        self._update_test_status(test_name, phase, pass_status if rc == 0 else TEST_FAIL_STATUS)
        return rc == 0

    ###########################################################################
    def _create_newcase(self, test_name):
    ###########################################################################
        print "Creating case for", test_name

        test_case, case_opts, grid, compset, machine, compiler, test_mods = acme_util.parse_test_name(test_name)
        scratch_dir = acme_util.get_machine_info("CESMSCRATCHROOT", machine=machine, project=self._project)

        test_dir = self._get_test_dir(test_name)

        create_newcase_cmd = "%s -silent -case %s -res %s -mach %s -compiler %s -compset %s -testname %s -project %s -nosavetiming -sharedlibroot %s" % \
                              (os.path.join(self._cime_root, "scripts", "create_newcase"),
                               test_dir, grid, machine, compiler, compset, test_case, self._project,
                               os.path.join(scratch_dir, "sharedlibroot.%s" % self._test_id))
        if (case_opts is not None):
            create_newcase_cmd += " -confopts _%s" % ("_".join(case_opts))
        if (test_mods is not None):
            test_mod_file = os.path.join(self._cime_root, "scripts", "Testing", "Testlistxml", "testmods_dirs", test_mods)
            if (not os.path.exists(test_mod_file)):
                self._log_output(test_name, "Missing testmod file '%s'" % test_mod_file)
                self._update_test_status(test_name, CREATE_NEWCASE_PHASE, TEST_FAIL_STATUS)
                return False
            create_newcase_cmd += " -user_mods_dir %s" % test_mod_file

        return self._run_phase_command(test_name, create_newcase_cmd, CREATE_NEWCASE_PHASE)

    ###########################################################################
    def _set_up_xml(self, test_name):
    ###########################################################################
        print "Setting up XML for", test_name

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
    def _setup_test(self, test_name):
    ###########################################################################
        print "Setting up test case", test_name

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
    def _process_namelists(self, test_name):
    ###########################################################################
        if (self._compare or self._generate):
            print "Processing namelists for test case", test_name

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
                self._tests_with_nl_problems.append(test_name)

        elif (self._generate):
            if (not os.path.isdir(baseline_dir)):
                os.makedirs(baseline_dir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH)

            if (os.path.isdir(baseline_casedocs)):
                shutil.rmtree(baseline_casedocs)
            shutil.copytree(casedoc_dir, baseline_casedocs)
            for item in glob.glob(os.path.join(test_dir, "user_nl*")):
                shutil.copy2(item, baseline_dir)

        # Always mark as passed unless we hit exception
        self._update_test_status(test_name, NAMELIST_PHASE, TEST_PASSED_STATUS)
        return True

    ###########################################################################
    def _build_test(self, test_name):
    ###########################################################################
        print "Building test case", test_name

        case_id = self._get_case_id(test_name)
        test_dir = self._get_test_dir(test_name)
        return self._run_phase_command(test_name, "./%s.test_build" % case_id, BUILD_PHASE, from_dir=test_dir)

    ###########################################################################
    def _run_test(self, test_name):
    ###########################################################################
        print "Running test case", test_name

        case_id = self._get_case_id(test_name)
        test_dir = self._get_test_dir(test_name)
        if (self._no_batch):
            return self._run_phase_command(test_name, "./%s.test" % case_id, RUN_PHASE, from_dir=test_dir, pass_status=TEST_PASSED_STATUS)
        else:
            return self._run_phase_command(test_name, "./%s.submit" % case_id, RUN_PHASE, from_dir=test_dir, pass_status=TEST_PENDING_STATUS)

    ###########################################################################
    def _test_status_phase(self, test_name):
    ###########################################################################
        str_to_write = ""
        made_it_to_phase = PHASES.index(self._get_test_phase(test_name))
        for phase in PHASES[0:made_it_to_phase+1]:
            str_to_write += "%s %s %s\n" % (self._get_test_status(test_name, phase), test_name, phase)

        if (not self._no_run and not self._is_broken(test_name)):
            # Ensure PEND state always gets added to TestStatus file if we are
            # about to run test
            str_to_write += "%s %s %s\n" % (TEST_PENDING_STATUS, test_name, RUN_PHASE)

        try:
            test_status_file = os.path.join(self._get_test_dir(test_name), TEST_STATUS_FILENAME)
            with open(test_status_file, "w") as fd:
                fd.write(str_to_write)
        except Exception as e:
            self._log_output(test_name, "VERY BAD! Could not make TestStatus file '%s': '%s'" % (test_status_file, str(e)))

    ###########################################################################
    def _run_catch_exceptions(self, phase, run):
    ###########################################################################
        for test_name in self._test_names:
            try:
                if (not self._is_broken(test_name)):
                    run(test_name)
            except Exception as e:
                exc_tb = sys.exc_info()[2]
                errput = "Test '%s' failed in phase '%s' with exception '%s'" % (test_name, phase, str(e))
                self._log_output(test_name, errput)
                warning("Caught exception: %s" % str(e))
                traceback.print_tb(exc_tb)
                self._update_test_status(test_name, phase, TEST_FAIL_STATUS)

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
        """
        # Tell user what will be run
        print "RUNNING TESTS:"
        for test_name in self._test_names:
            print " ", test_name

        # TODO - do in parallel
        # TODO - unit tests where possible
        # TODO - documentation

        # Create cases
        self._run_catch_exceptions(CREATE_NEWCASE_PHASE, self._create_newcase)

        # Set up XML
        self._run_catch_exceptions(XML_PHASE, self._set_up_xml)

        # Setup
        self._run_catch_exceptions(SETUP_PHASE, self._setup_test)

        # Namelists
        self._run_catch_exceptions(NAMELIST_PHASE, self._process_namelists)

        # Build
        if (not self._no_build):
            self._run_catch_exceptions(BUILD_PHASE, self._build_test)

        # TODO: The run scripts heavily use the TestStatus file. So we write out
        # the phases we have taken care of and then let the run scrips go from there
        # Eventually, it would be nice to have TestStatus management encapsulated
        # into a single place.

        # Make status files
        for test_name in self._test_names:
            self._test_status_phase(test_name)

        # Run
        if (not self._no_run):
            self._run_catch_exceptions(RUN_PHASE, self._run_test)

        # Setup cs files
        self._setup_cs_files()

        # If we failed VERY early on in the run phase, it's possible that
        # the CIME scripts never got a chance to set the state.
        if (not self._no_run):
            for test_name in self._test_names:
                if (self._get_test_phase(test_name) == RUN_PHASE and
                    self._is_broken(test_name)):
                    try:
                        test_status_file = os.path.join(self._get_test_dir(test_name), TEST_STATUS_FILENAME)
                        statuses = wait_for_tests.parse_test_status_file(test_status_file)[0]
                        if (RUN_PHASE not in statuses):
                            self._test_status_phase(test_name)
                        elif (statuses[RUN_PHASE] in [TEST_PASSED_STATUS, TEST_PENDING_STATUS]):
                            self._log_output(test_name,
                                             "VERY BAD! How was infrastructure able to log a TestState but not change it to FAIL?")
                    except Exception as e:
                        self._log_output(test_name, "VERY BAD! Could not read TestStatus file '%s': '%s'" % (test_status_file, str(e)))

        # Return True if all tests passed
        print "At create_test close, state is:"
        rv = True
        for idx, test_name in enumerate(self._test_names):
            phase, status = self._test_states[idx]
            if (status not in [TEST_PASSED_STATUS, TEST_PENDING_STATUS]):
                print "%s %s (phase %s)" % (status, test_name, phase)
                rv = False
            elif (test_name in self._tests_with_nl_problems):
                print "%s %s (but otherwise OK)" % (NAMELIST_FAIL_STATUS, test_name)
                rv = False
            else:
                print status, test_name

        return rv
