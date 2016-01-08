"""
Encapsulate the importing of python utils
"""

import sys, os

_LIB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils", "python")
sys.path.append(_LIB_DIR)

print _LIB_DIR
