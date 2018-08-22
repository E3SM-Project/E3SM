"""
Implementation of the cs.status script, which prints the status of all
of the tests in one or more test suites
"""

from __future__ import print_function
from CIME.XML.standard_module_setup import *
from CIME.test_status import TestStatus
import os
import sys

def cs_status(test_paths, summary=False, out=sys.stdout):
    """Print the test statuses of all tests in test_paths. The default
    is to print to stdout, but this can be overridden with the 'out'
    argument.
    """
    test_id_output = {}
    for test_path in test_paths:
        test_dir=os.path.dirname(test_path)
        ts = TestStatus(test_dir=test_dir)
        test_id = os.path.basename(test_dir).split(".")[-1]
        test_name = ts.get_name()
        status = ts.get_overall_test_status()
        if not summary:
            output = "  %s (Overall: %s) details:\n" % (test_name, status)
            output += ts.phase_statuses_dump(prefix="    ")
        else:
            output = "  %s %s\n" % (status, test_name)

        if test_id in test_id_output:
            test_id_output[test_id] += output
        else:
            test_id_output[test_id] = output

    for test_id in sorted(test_id_output):
        print(test_id, file=out)
        print(test_id_output[test_id], file=out)
        print(' ', file=out)
