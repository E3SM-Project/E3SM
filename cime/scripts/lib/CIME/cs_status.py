"""
Implementation of the cs.status script, which prints the status of all
of the tests in one or more test suites
"""

from __future__ import print_function
from CIME.XML.standard_module_setup import *
from CIME.XML.expected_fails_file import ExpectedFailsFile
from CIME.test_status import TestStatus
import os
import sys
from collections import defaultdict

def cs_status(test_paths, summary=False, fails_only=False,
              count_fails_phase_list=None,
              expected_fails_filepath=None,
              out=sys.stdout):
    """Print the test statuses of all tests in test_paths. The default
    is to print to stdout, but this can be overridden with the 'out'
    argument.

    If summary is True, then only the overall status of each test is printed

    If fails_only is True, then only test failures are printed (this
    includes PENDs as well as FAILs).

    If count_fails_phase_list is provided, it should be a list of phases
    (from the phases given by test_status.ALL_PHASES). For each phase in
    this list: do not give line-by-line output; instead, just report the
    total number of tests that have not PASSed this phase (this includes
    PENDs and FAILs). (This is typically used with the fails_only
    option, but it can also be used without that option.)

    If expected_fails_filepath is provided, it should be a string giving
    the full path to a file listing expected failures for this test
    suite. Expected failures are then labeled as such in the output.
    """
    expect(not (summary and fails_only),
           "Cannot have both summary and fails_only")
    expect(not (summary and count_fails_phase_list),
           "Cannot have both summary and count_fails_phase_list")
    if count_fails_phase_list is None:
        count_fails_phase_list = []
    non_pass_counts = dict.fromkeys(count_fails_phase_list, 0)
    xfails = _get_xfails(expected_fails_filepath)
    test_id_output = defaultdict(str)
    test_id_counts = defaultdict(int)
    for test_path in test_paths:
        test_dir=os.path.dirname(test_path)
        ts = TestStatus(test_dir=test_dir)
        test_id = os.path.basename(test_dir).split(".")[-1]
        if summary:
            output = _overall_output(ts, "  {status} {test_name}\n")
        else:
            if fails_only:
                output = ''
            else:
                output = _overall_output(ts, "  {test_name} (Overall: {status}) details:\n")
            output += ts.phase_statuses_dump(prefix="    ",
                                             skip_passes=fails_only,
                                             skip_phase_list=count_fails_phase_list,
                                             xfails=xfails.get(ts.get_name()))
            if count_fails_phase_list:
                ts.increment_non_pass_counts(non_pass_counts)

        test_id_output[test_id] += output
        test_id_counts[test_id] += 1

    for test_id in sorted(test_id_output):
        count = test_id_counts[test_id]
        print("{}: {} test{}".format(test_id, count, 's' if count > 1 else ''), file=out)
        print(test_id_output[test_id], file=out)
        print(' ', file=out)

    if count_fails_phase_list:
        print(72*'=', file=out)
        print('Non-PASS results for select phases:', file=out)
        for phase in count_fails_phase_list:
            print('{} non-passes: {}'.format(phase, non_pass_counts[phase]), file=out)

def _get_xfails(expected_fails_filepath):
    """Returns a dictionary of ExpectedFails objects, where the keys are test names

    expected_fails_filepath should be either a string giving the path to
    the file containing expected failures, or None. If None, then this
    returns an empty dictionary (as if expected_fails_filepath were
    pointing to a file with no expected failures listed).
    """
    if expected_fails_filepath is not None:
        expected_fails_file = ExpectedFailsFile(expected_fails_filepath)
        xfails = expected_fails_file.get_expected_fails()
    else:
        xfails = {}
    return xfails

def _overall_output(ts, format_str):
    """Returns a string giving the overall test status

    Args:
    ts: TestStatus object
    format_str (string): string giving the format of the output; must
        contain place-holders for status and test_name
    """
    test_name = ts.get_name()
    status = ts.get_overall_test_status()
    return format_str.format(status=status, test_name=test_name)
