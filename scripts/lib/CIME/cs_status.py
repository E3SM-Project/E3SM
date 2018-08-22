"""
Implementation of the cs.status script, which prints the status of all
of the tests in one or more test suites
"""

from __future__ import print_function
from CIME.XML.standard_module_setup import *
from CIME.test_status import TestStatus
import os
import sys

def cs_status(test_paths, summary=False, fails_only=False, out=sys.stdout):
    """Print the test statuses of all tests in test_paths. The default
    is to print to stdout, but this can be overridden with the 'out'
    argument.

    If summary is True, then only the overall status of each test is printed

    If fails_only is True, then only test failures are printed (this
    includes PENDs as well as FAILs).
    """
    expect(not (summary and fails_only), "Cannot have both summary and fails_only")
    test_id_output = {}
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
            output += ts.phase_statuses_dump(prefix="    ", skip_passes=fails_only)

        if test_id in test_id_output:
            test_id_output[test_id] += output
        else:
            test_id_output[test_id] = output

    for test_id in sorted(test_id_output):
        print(test_id, file=out)
        print(test_id_output[test_id], file=out)
        print(' ', file=out)

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
