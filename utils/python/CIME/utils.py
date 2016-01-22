"""
Common functions used by cime python scripts
"""

import sys
import os

def expect(condition, error_msg):
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.

    >>> expect(True, "error1")
    >>> expect(False, "error2")
    Traceback (most recent call last):
        ...
    SystemExit: FAIL: error2
    """
    if (not condition):
        raise SystemExit("FAIL: %s" % error_msg)
