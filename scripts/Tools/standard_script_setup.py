"""
Encapsulate the importing of python utils and logging setup, things
that every script should do.
"""
# pylint: disable=unused-import

import sys, os
import __main__ as main
_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..")
_LIB_DIR = os.path.join(_CIMEROOT, "utils", "python")
sys.path.append(_LIB_DIR)

# Important: Allows external tools to link up with CIME
os.environ["CIMEROOT"] = _CIMEROOT

import CIME.utils
CIME.utils.check_minimum_python_version(2, 7)
CIME.utils.stop_buffering_output()
import logging, doctest, argparse
