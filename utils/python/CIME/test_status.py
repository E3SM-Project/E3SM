"""
Contains the crucial TestStatus class which manages phase-state of a test
case and ensure that this state is represented by the TestStatus file in
the case.

TestStatus objects are only modifiable via the set_status method and this
is only allowed if the object is being accessed within the context of a
context manager. Example:

    with TestStatus(test_dir=caseroot) as ts:
        ts.set_status(RUN_PHASE, TEST_PASS_STATUS)

This file also contains all of the hardcoded phase information which includes
the phase names, phase orders, potential phase states, and which phases are
required (core phases).

Additional important design decisions:
1) In order to ensure that incomplete tests are always left in a PEND
   state, updating a core phase to a PASS state will automatically set the next
   core state to PEND.
2) If the user repeats a core state, that invalidates all subsequent state. For
   example, if a user rebuilds their case, then any of the post-run states like the
   RUN state are no longer valid.
"""

from CIME.XML.standard_module_setup import *

from collections import OrderedDict

TEST_STATUS_FILENAME = "TestStatus"

# The statuses that a phase can be in
TEST_PEND_STATUS = "PEND"
TEST_RUN_STATUS = "RUN"
TEST_PASS_STATUS = "PASS"
TEST_FAIL_STATUS = "FAIL"

ALL_PHASE_STATUSES = [TEST_PEND_STATUS, TEST_RUN_STATUS, TEST_PASS_STATUS, TEST_FAIL_STATUS]

# Special statuses that the overall test can be in
TEST_DIFF_STATUS     = "DIFF"   # Implies a failure in the BASELINE phase
NAMELIST_FAIL_STATUS = "NLFAIL" # Implies a failure in the NLCOMP phase

# Special strings that can appear in comments, indicating particular types of failures
TEST_NO_BASELINES_COMMENT = "BFAIL" # Implies baseline directory is missing in the baseline comparison phase

# The valid phases
CREATE_NEWCASE_PHASE  = "CREATE_NEWCASE"
XML_PHASE             = "XML"
SETUP_PHASE           = "SETUP"
NAMELIST_PHASE        = "NLCOMP"
SHAREDLIB_BUILD_PHASE = "SHAREDLIB_BUILD"
MODEL_BUILD_PHASE     = "MODEL_BUILD"
RUN_PHASE             = "RUN"
THROUGHPUT_PHASE      = "TPUTCOMP"
MEMCOMP_PHASE         = "MEMCOMP"
MEMLEAK_PHASE         = "MEMLEAK"
COMPARE_PHASE         = "COMPARE" # This is one special, real phase will be COMPARE_$WHAT, this is for internal test comparisons, there could be multiple variations of this phase in one test
BASELINE_PHASE        = "BASELINE"
GENERATE_PHASE        = "GENERATE"

ALL_PHASES = [CREATE_NEWCASE_PHASE,
              XML_PHASE,
              SETUP_PHASE,
              NAMELIST_PHASE,
              SHAREDLIB_BUILD_PHASE,
              MODEL_BUILD_PHASE,
              RUN_PHASE,
              COMPARE_PHASE,
              BASELINE_PHASE,
              THROUGHPUT_PHASE,
              MEMCOMP_PHASE,
              MEMLEAK_PHASE,
              GENERATE_PHASE]

# These are mandatory phases that a test must go through
CORE_PHASES = [CREATE_NEWCASE_PHASE,
               XML_PHASE,
               SETUP_PHASE,
               SHAREDLIB_BUILD_PHASE,
               MODEL_BUILD_PHASE,
               RUN_PHASE]

def _test_helper1(file_contents):
    ts = TestStatus(test_dir="/", test_name="ERS.foo.A")
    ts._parse_test_status(file_contents) # pylint: disable=protected-access
    return ts._phase_statuses # pylint: disable=protected-access

def _test_helper2(file_contents, wait_for_run=False, check_throughput=False, check_memory=False, ignore_namelists=False):
    ts = TestStatus(test_dir="/", test_name="ERS.foo.A")
    ts._parse_test_status(file_contents) # pylint: disable=protected-access
    return ts.get_overall_test_status(wait_for_run=wait_for_run,
                                      check_throughput=check_throughput,
                                      check_memory=check_memory,
                                      ignore_namelists=ignore_namelists)

class TestStatus(object):

    def __init__(self, test_dir=None, test_name=None, no_io=False):
        """
        Create a TestStatus object

        If test_dir is not specified, it is set to the current working directory

        no_io is intended only for testing, and should be kept False in
        production code
        """
        test_dir = os.getcwd() if test_dir is None else test_dir
        self._filename = os.path.join(test_dir, TEST_STATUS_FILENAME)
        self._phase_statuses = OrderedDict() # {name -> (status, comments)}
        self._test_name = test_name
        self._ok_to_modify = False
        self._no_io = no_io

        if os.path.exists(self._filename):
            self._parse_test_status_file()
        else:
            expect(test_name is not None, "Must provide test_name if TestStatus file doesn't exist")

    def __enter__(self):
        self._ok_to_modify = True
        return self

    def __exit__(self, *_):
        self._ok_to_modify = False
        self.flush()

    def __iter__(self):
        for phase, data in self._phase_statuses.iteritems():
            yield phase, data[0]

    def get_name(self):
        return self._test_name

    def set_status(self, phase, status, comments=""):
        """
        Update the status of this test by changing the status of given phase to the
        given status.

        >>> with TestStatus(test_dir="/", test_name="ERS.foo.A", no_io=True) as ts:
        ...     ts.set_status(CREATE_NEWCASE_PHASE, "PASS")
        ...     ts.set_status(XML_PHASE, "PASS")
        ...     ts.set_status(SETUP_PHASE, "FAIL")
        ...     ts.set_status(SETUP_PHASE, "PASS")
        ...     ts.set_status("%s_base_rest" % COMPARE_PHASE, "FAIL")
        ...     ts.set_status(SHAREDLIB_BUILD_PHASE, "PASS", comments='Time=42')
        >>> ts._phase_statuses
        OrderedDict([('CREATE_NEWCASE', ('PASS', '')), ('XML', ('PASS', '')), ('SETUP', ('PASS', '')), ('SHAREDLIB_BUILD', ('PASS', 'Time=42')), ('COMPARE_base_rest', ('FAIL', '')), ('MODEL_BUILD', ('PEND', ''))])

        >>> with TestStatus(test_dir="/", test_name="ERS.foo.A", no_io=True) as ts:
        ...     ts.set_status(CREATE_NEWCASE_PHASE, "PASS")
        ...     ts.set_status(XML_PHASE, "PASS")
        ...     ts.set_status(SETUP_PHASE, "FAIL")
        ...     ts.set_status(SETUP_PHASE, "PASS")
        ...     ts.set_status(BASELINE_PHASE, "PASS")
        ...     ts.set_status("%s_base_rest" % COMPARE_PHASE, "FAIL")
        ...     ts.set_status(SHAREDLIB_BUILD_PHASE, "PASS", comments='Time=42')
        ...     ts.set_status(SETUP_PHASE, "PASS")
        >>> ts._phase_statuses
        OrderedDict([('CREATE_NEWCASE', ('PASS', '')), ('XML', ('PASS', '')), ('SETUP', ('PASS', '')), ('SHAREDLIB_BUILD', ('PEND', ''))])

        >>> with TestStatus(test_dir="/", test_name="ERS.foo.A", no_io=True) as ts:
        ...     ts.set_status(CREATE_NEWCASE_PHASE, "FAIL")
        >>> ts._phase_statuses
        OrderedDict([('CREATE_NEWCASE', ('FAIL', ''))])
        """
        expect(self._ok_to_modify, "TestStatus not in a modifiable state, use 'with' syntax")
        expect(phase in ALL_PHASES or phase.startswith(COMPARE_PHASE),
               "Invalid phase '%s'" % phase)
        expect(status in ALL_PHASE_STATUSES, "Invalid status '%s'" % status)

        if phase in CORE_PHASES and phase != CORE_PHASES[0]:
            previous_core_phase = CORE_PHASES[CORE_PHASES.index(phase)-1]
            #TODO: enamble check below
            #expect(previous_core_phase in self._phase_statuses, "Core phase '%s' was skipped" % previous_core_phase)

            if previous_core_phase in self._phase_statuses:
                expect(self._phase_statuses[previous_core_phase][0] == TEST_PASS_STATUS,
                       "Cannot move past core phase '%s', it didn't pass" % previous_core_phase)

        reran_phase = (phase in self._phase_statuses and self._phase_statuses[phase][0] != TEST_PEND_STATUS and phase in CORE_PHASES)
        if reran_phase:
            # All subsequent phases are invalidated
            phase_idx = ALL_PHASES.index(phase)
            for subsequent_phase in ALL_PHASES[phase_idx+1:]:
                if subsequent_phase in self._phase_statuses:
                    del self._phase_statuses[subsequent_phase]
                if subsequent_phase.startswith(COMPARE_PHASE):
                    for stored_phase in self._phase_statuses.keys():
                        if stored_phase.startswith(COMPARE_PHASE):
                            del self._phase_statuses[stored_phase]

        self._phase_statuses[phase] = (status, comments) # Can overwrite old phase info

        if status == TEST_PASS_STATUS and phase in CORE_PHASES and phase != CORE_PHASES[-1]:
            next_core_phase = CORE_PHASES[CORE_PHASES.index(phase)+1]
            self._phase_statuses[next_core_phase] = (TEST_PEND_STATUS, "")

    def get_status(self, phase):
        return self._phase_statuses[phase][0] if phase in self._phase_statuses else None

    def get_comment(self, phase):
        return self._phase_statuses[phase][1] if phase in self._phase_statuses else None

    def phase_statuses_dump(self, fd, prefix=''):
        """
        Args:
            fd: file open for writing
            prefix: string printed at the start of each line
        """
        if self._phase_statuses:
            for phase, data in self._phase_statuses.iteritems():
                status, comments = data
                if not comments:
                    fd.write("%s%s %s %s\n" % (prefix, status, self._test_name, phase))
                else:
                    fd.write("%s%s %s %s %s\n" % (prefix, status, self._test_name, phase, comments))

    def flush(self):
        if self._phase_statuses and not self._no_io:
            with open(self._filename, "w") as fd:
                self.phase_statuses_dump(fd)

    def _parse_test_status(self, file_contents):
        """
        >>> contents = '''
        ... PASS ERS.foo.A CREATE_NEWCASE
        ... PASS ERS.foo.A XML
        ... FAIL ERS.foo.A SETUP
        ... PASS ERS.foo.A COMPARE_base_rest
        ... PASS ERS.foo.A SHAREDLIB_BUILD Time=42
        ... '''
        >>> _test_helper1(contents)
        OrderedDict([('CREATE_NEWCASE', ('PASS', '')), ('XML', ('PASS', '')), ('SETUP', ('FAIL', '')), ('COMPARE_base_rest', ('PASS', '')), ('SHAREDLIB_BUILD', ('PASS', 'Time=42'))])
        """
        for line in file_contents.splitlines():
            line = line.strip()
            tokens = line.split()
            if line == "":
                pass # skip blank lines
            elif len(tokens) >= 3:
                status, curr_test_name, phase = tokens[:3]
                if (self._test_name is None):
                    self._test_name = curr_test_name
                else:
                    expect(self._test_name == curr_test_name,
                           "inconsistent test name in parse_test_status: '%s' != '%s'" % (self._test_name, curr_test_name))

                expect(status in ALL_PHASE_STATUSES,
                       "Unexpected status '%s' in parse_test_status for test '%s'" % (status, self._test_name))
                expect(phase in ALL_PHASES or phase.startswith(COMPARE_PHASE),
                       "phase '%s' not expected in parse_test_status for test '%s'" % (phase, self._test_name))
                expect(phase not in self._phase_statuses,
                       "Should not have seen multiple instances of phase '%s' for test '%s'" % (phase, self._test_name))

                self._phase_statuses[phase] = (status, " ".join(tokens[3:]))
            else:
                logging.warning("In TestStatus file for test '%s', line '%s' not in expected format" % (self._test_name, line))

    def _parse_test_status_file(self):
        with open(self._filename, "r") as fd:
            self._parse_test_status(fd.read())

    def get_overall_test_status(self, wait_for_run=False, check_throughput=False, check_memory=False, ignore_namelists=False, ignore_memleak=False):
        r"""
        Given the current phases and statuses, produce a single results for this test. Preference
        is given to PEND since we don't want to stop waiting for a test
        that hasn't finished. Namelist diffs are given the lowest precedence.

        >>> _test_helper2('PASS ERS.foo.A RUN')
        'PASS'
        >>> _test_helper2('PASS ERS.foo.A SHAREDLIB_BUILD\nPEND ERS.foo.A RUN')
        'PEND'
        >>> _test_helper2('FAIL ERS.foo.A MODEL_BUILD\nPEND ERS.foo.A RUN')
        'PEND'
        >>> _test_helper2('PASS ERS.foo.A MODEL_BUILD\nPASS ERS.foo.A RUN')
        'PASS'
        >>> _test_helper2('PASS ERS.foo.A RUN\nFAIL ERS.foo.A TPUTCOMP')
        'PASS'
        >>> _test_helper2('PASS ERS.foo.A RUN\nFAIL ERS.foo.A TPUTCOMP', check_throughput=True)
        'FAIL'
        >>> _test_helper2('PASS ERS.foo.A RUN\nFAIL ERS.foo.A NLCOMP')
        'NLFAIL'
        >>> _test_helper2('PASS ERS.foo.A RUN\nFAIL ERS.foo.A MEMCOMP')
        'PASS'
        >>> _test_helper2('PASS ERS.foo.A RUN\nFAIL ERS.foo.A NLCOMP', ignore_namelists=True)
        'PASS'
        >>> _test_helper2('PASS ERS.foo.A COMPARE_1\nFAIL ERS.foo.A NLCOMP\nFAIL ERS.foo.A COMPARE_2\nPASS ERS.foo.A RUN')
        'FAIL'
        >>> _test_helper2('FAIL ERS.foo.A BASELINE\nFAIL ERS.foo.A NLCOMP\nPASS ERS.foo.A COMPARE_2\nPASS ERS.foo.A RUN')
        'DIFF'
        >>> _test_helper2('FAIL ERS.foo.A BASELINE\nFAIL ERS.foo.A NLCOMP\nFAIL ERS.foo.A COMPARE_2\nPASS ERS.foo.A RUN')
        'FAIL'
        >>> _test_helper2('PASS ERS.foo.A MODEL_BUILD')
        'PASS'
        >>> _test_helper2('PASS ERS.foo.A MODEL_BUILD', wait_for_run=True)
        'PEND'
        >>> _test_helper2('FAIL ERS.foo.A MODEL_BUILD', wait_for_run=True)
        'FAIL'
        >>> _test_helper2('PASS ERS.foo.A MODEL_BUILD\nPEND ERS.foo.A RUN', wait_for_run=True)
        'PEND'
        >>> _test_helper2('PASS ERS.foo.A MODEL_BUILD\nFAIL ERS.foo.A RUN', wait_for_run=True)
        'FAIL'
        >>> _test_helper2('PASS ERS.foo.A MODEL_BUILD\nPASS ERS.foo.A RUN', wait_for_run=True)
        'PASS'
        """
        rv = TEST_PASS_STATUS
        run_phase_found = False
        for phase, data in self._phase_statuses.iteritems():
            status = data[0]
            if phase == RUN_PHASE:
                run_phase_found = True

            if (status == TEST_PEND_STATUS):
                return status

            elif (status == TEST_FAIL_STATUS):
                if ( (not check_throughput and phase == THROUGHPUT_PHASE) or
                     (not check_memory and phase == MEMCOMP_PHASE) or
                     (ignore_namelists and phase == NAMELIST_PHASE) or
                     (ignore_memleak and phase == MEMLEAK_PHASE) ):
                    continue

                if (phase == NAMELIST_PHASE):
                    if (rv == TEST_PASS_STATUS):
                        rv = NAMELIST_FAIL_STATUS

                elif (rv in [NAMELIST_FAIL_STATUS, TEST_PASS_STATUS] and phase == BASELINE_PHASE):
                    rv = TEST_DIFF_STATUS

                else:
                    rv = TEST_FAIL_STATUS

        # The test did not fail but the RUN phase was not found, so if the user requested
        # that we wait for the RUN phase, then the test must still be considered pending.
        if rv != TEST_FAIL_STATUS and not run_phase_found and wait_for_run:
            rv = TEST_PEND_STATUS

        return rv
