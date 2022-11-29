
import sys, os

# Add CIME libs to sys path
_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","cime")
_LIB_DIR = os.path.join(_CIMEROOT)
sys.path.append(_LIB_DIR)

from CIME.XML.machines import Machines

###############################################################################
def query_cime(machine, param):
###############################################################################
    mach_obj = Machines(machine=machine)
    return mach_obj.get_value(param)
