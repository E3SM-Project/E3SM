"""
Encapsulate the importing of python utils and logging setup, things
that every script should do.
"""

import sys, os, logging, doctest, argparse, logging.config
import __main__ as main
_CIMEROOT = os.environ.get("CIMEROOT")
if(_CIMEROOT is None):
    _CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..")
    os.environ["CIMEROOT"] = _CIMEROOT
_LIB_DIR = os.path.join(_CIMEROOT, "utils", "python")
sys.path.append(_LIB_DIR)


import CIME.utils
CIME.utils.check_minimum_python_version(2, 7)
CIME.utils.stop_buffering_output()
# set up logging to file
logging.basicConfig(level=logging.INFO)
