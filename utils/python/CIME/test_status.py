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
MEMORY_PHASE          = "MEMCOMP"
COMPARE_PHASE         = "COMPARE"

ALL_PHASES = [INITIAL_PHASE,
              CREATE_NEWCASE_PHASE,
              XML_PHASE,
              SETUP_PHASE,
              NAMELIST_PHASE,
              SHAREDLIB_BUILD_PHASE,
              MODEL_BUILD_PHASE,
              RUN_PHASE,
              THROUGHPUT_PHASE,
              MEMORY_PHASE,
              COMPARE_PHASE]

MULTI_PHASES = [COMPARE_PHASE]

class TestStatus(object):

    def __init__(self, test_dir=os.getcwd()):
        self._filename = os.path.join(test_dir, TEST_STATUS_FILENAME)
        self._phase_statuses = OrderedDict() # {name -> (status, comments)}
        self._test_name = None

        if os.path.exists(self._filename):
            self._parse_test_status_file()

    def _parse_test_status(self, file_contents):
        """
        >>> ts = TestStatus()
        >>> contents = '''
        ... CREATE_NEWCASE ERS.foo.A PASS
        ... XML_PHASE      ERS.foo.A PASS
        ... SETUP_PHASE    ERS.foo.A FAIL
        ... '''
        >>> ts._parse_test_status(contents)
        >>> ts._phase_statuses
        '{stuff}'
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
                    expect(self._test_name == curr_test_name, "inconsistent test name in parse_test_status: '%s' != '%s'"%(self._test_name, curr_test_name))

                expect(status in ALL_PHASE_STATUSES,
                       "Unexpected status '%s' in parse_test_status for test '%s'" % (status, self._test_name))
                expect(phase in ALL_PHASES,
                       "phase '%s' not expected in parse_test_status for test '%s'" % (phase, self._test_name))

                if (phase in rv):
                    # Phase names don't matter here, just need something unique
                    rv[phase] = reduce_stati({"%s_" % phase : status, phase : rv[phase]})
                else:
                    rv[phase] = status
            else:
                logging.warning("In TestStatus file for test '%s', line '%s' not in expected format" % (self._test_name, line))

    def _parse_test_status_file(self):
        with open(self._filename, "r") as fd:
            self._parse_test_status(fd.read())


def reduce_stati(stati, wait_for_run=False, check_throughput=False, check_memory=False, ignore_namelists=False):
    """
    Given a collection of stati for a test, produce a single result. Preference
    is given to unfinished stati since we don't want to stop waiting for a test
    that hasn't finished. Namelist diffs are given the lowest precedence.
    """
    rv = TEST_PASS_STATUS
    run_phase_found = False
    for phase, status in stati.iteritems():
        if phase == RUN_PHASE:
            run_phase_found = True

        if (status == TEST_PENDING_STATUS):
            return status

        elif (status != TEST_PASS_STATUS):
            if ( (not check_throughput and THROUGHPUT_TEST_STR in phase) or
                 (not check_memory and MEMORY_TEST_STR in phase) or
                 (ignore_namelists and phase == NAMELIST_PHASE) ):
                continue

            if (status == NAMELIST_FAIL_STATUS):
                if (rv == TEST_PASS_STATUS):
                    rv = NAMELIST_FAIL_STATUS

            elif (rv in [NAMELIST_FAIL_STATUS, TEST_PASS_STATUS] and phase == HIST_COMPARE_PHASE):
                rv = TEST_DIFF_STATUS

            else:
                rv = status

    # The test did not fail but the RUN phase was not found, so if the user requested
    # that we wait for the RUN phase, then the test must still be considered pending.
    if rv != TEST_FAIL_STATUS and not run_phase_found and wait_for_run:
        rv = TEST_PENDING_STATUS

    return rv


def interpret_status(file_contents, wait_for_run=False, check_throughput=False, check_memory=False, ignore_namelists=False):
    r"""
    >>> interpret_status('PASS testname RUN')
    ('testname', 'PASS')
    >>> interpret_status('PASS testname SHAREDLIB_BUILD\nPEND testname RUN')
    ('testname', 'PEND')
    >>> interpret_status('FAIL testname MODEL_BUILD\nPEND testname RUN')
    ('testname', 'PEND')
    >>> interpret_status('PASS testname MODEL_BUILD\nPASS testname RUN')
    ('testname', 'PASS')
    >>> interpret_status('PASS testname RUN\nFAIL testname tputcomp')
    ('testname', 'PASS')
    >>> interpret_status('PASS testname RUN\nFAIL testname tputcomp', check_throughput=True)
    ('testname', 'FAIL')
    >>> interpret_status('PASS testname RUN\nNLFAIL testname nlcomp')
    ('testname', 'NLFAIL')
    >>> interpret_status('PASS testname RUN\nFAIL testname memleak')
    ('testname', 'FAIL')
    >>> interpret_status('PASS testname RUN\nNLFAIL testname nlcomp', ignore_namelists=True)
    ('testname', 'PASS')
    >>> interpret_status('PASS testname compare\nNLFAIL testname nlcomp\nFAIL testname compare\nPASS testname RUN')
    ('testname', 'DIFF')
    >>> interpret_status('PASS testname MODEL_BUILD')
    ('testname', 'PASS')
    >>> interpret_status('PASS testname MODEL_BUILD', wait_for_run=True)
    ('testname', 'PEND')
    """
    statuses, test_name = parse_test_status(file_contents)
    reduced_status = reduce_stati(statuses, wait_for_run, check_throughput, check_memory, ignore_namelists)

    return test_name, reduced_status

def interpret_status_file(file_name, wait_for_run=False, check_throughput=False, check_memory=False, ignore_namelists=False):
    with open(file_name, "r") as fd:
        return interpret_status(fd.read(), wait_for_run, check_throughput, check_memory, ignore_namelists)
