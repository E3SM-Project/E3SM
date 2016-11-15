"""
Functions for managing the TestStatus file
"""

from CIME.XML.standard_module_setup import *

from collections import OrderedDict

TEST_STATUS_FILENAME = "TestStatus"

# The statuses that a phase can be in
TEST_PENDING_STATUS  = "PEND"
TEST_PASS_STATUS     = "PASS"
TEST_FAIL_STATUS     = "FAIL"

ALL_PHASE_STATUSES = [TEST_PENDING_STATUS, TEST_PASS_STATUS, TEST_FAIL_STATUS]

# Special statuses that the overall test can be in
TEST_DIFF_STATUS     = "DIFF"   # Implies a failure in one of the COMPARE phases
NAMELIST_FAIL_STATUS = "NLFAIL" # Implies a failure in the NLCOMP phase

# The valid phases
INITIAL_PHASE         = "INIT"
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
COMPARE_PHASE         = "COMPARE" # This is one special, real phase will be COMPARE_$WHAT
GENERATE_PHASE        = "GENERATE"

ALL_PHASES = [INITIAL_PHASE,
              CREATE_NEWCASE_PHASE,
              XML_PHASE,
              SETUP_PHASE,
              NAMELIST_PHASE,
              SHAREDLIB_BUILD_PHASE,
              MODEL_BUILD_PHASE,
              RUN_PHASE,
              THROUGHPUT_PHASE,
              MEMCOMP_PHASE,
              MEMLEAK_PHASE,
              GENERATE_PHASE]

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

    def __init__(self, test_dir=None, test_name=None):
        """
        Create a TestStatus object

        If test_dir is not specified, it is set to the current working directory
        """
        test_dir = os.getcwd() if test_dir is None else test_dir
        self._filename = os.path.join(test_dir, TEST_STATUS_FILENAME)
        self._phase_statuses = OrderedDict() # {name -> (status, comments)}
        self._test_name = test_name
        self._ok_to_modify = False

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
        >>> with TestStatus(test_dir="/", test_name="ERS.foo.A") as ts:
        ...     ts.set_status(CREATE_NEWCASE_PHASE, "PASS")
        ...     ts.set_status(XML_PHASE, "PASS")
        ...     ts.set_status(SETUP_PHASE, "PASS")
        ...     ts.set_status(SETUP_PHASE, "FAIL")
        ...     ts.set_status("%s_baseline" % COMPARE_PHASE, "FAIL")
        ...     ts.set_status(SHAREDLIB_BUILD_PHASE, "PASS", comments='Time=42')
        ...     result = OrderedDict(ts._phase_statuses)
        ...     ts._phase_statuses = None
        >>> result
        OrderedDict([('CREATE_NEWCASE', ('PASS', '')), ('XML', ('PASS', '')), ('SETUP', ('FAIL', '')), ('COMPARE_baseline', ('FAIL', '')), ('SHAREDLIB_BUILD', ('PASS', 'Time=42'))])
        """
        expect(self._ok_to_modify, "TestStatus not in a modifiable state, use 'with' syntax")
        expect(phase in ALL_PHASES or phase.startswith(COMPARE_PHASE),
               "Invalid phase '%s'" % phase)
        expect(status in ALL_PHASE_STATUSES, "Invalid status '%s'" % status)

        self._phase_statuses[phase] = (status, comments) # Can overwrite old phase info

    def get_status(self, phase):
        return self._phase_statuses[phase][0] if phase in self._phase_statuses else None

    def get_comment(self, phase):
        return self._phase_statuses[phase][1] if phase in self._phase_statuses else None

    def flush(self):
        if self._phase_statuses:
            with open(self._filename, "w") as fd:
                for phase, data in self._phase_statuses.iteritems():
                    status, comments = data
                    if not comments:
                        fd.write("%s %s %s\n" % (status, self._test_name, phase))
                    else:
                        fd.write("%s %s %s %s\n" % (status, self._test_name, phase, comments))

    def _parse_test_status(self, file_contents):
        """
        >>> contents = '''
        ... PASS ERS.foo.A CREATE_NEWCASE
        ... PASS ERS.foo.A XML
        ... FAIL ERS.foo.A SETUP
        ... PASS ERS.foo.A COMPARE_baseline
        ... PASS ERS.foo.A SHAREDLIB_BUILD Time=42
        ... '''
        >>> _test_helper1(contents)
        OrderedDict([('CREATE_NEWCASE', ('PASS', '')), ('XML', ('PASS', '')), ('SETUP', ('FAIL', '')), ('COMPARE_baseline', ('PASS', '')), ('SHAREDLIB_BUILD', ('PASS', 'Time=42'))])
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
        'DIFF'
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

            if (status == TEST_PENDING_STATUS):
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

                elif (rv in [NAMELIST_FAIL_STATUS, TEST_PASS_STATUS] and phase.startswith(COMPARE_PHASE)):
                    rv = TEST_DIFF_STATUS

                else:
                    rv = TEST_FAIL_STATUS

        # The test did not fail but the RUN phase was not found, so if the user requested
        # that we wait for the RUN phase, then the test must still be considered pending.
        if rv != TEST_FAIL_STATUS and not run_phase_found and wait_for_run:
            rv = TEST_PENDING_STATUS

        return rv
