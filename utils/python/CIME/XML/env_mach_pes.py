"""
Interface to the env_mach_pes.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvMachPes(EnvBase):

    def __init__(self, case_root=os.getcwd(), infile="env_mach_pes.xml"):
        """
        initialize an object interface to file env_mach_pes.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)

    def get_value(self, vid, attribute={}, resolved=True, subgroup=None):
        value = EnvBase.get_value(self, vid, attribute, resolved, subgroup)
        if "NTASKS" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        if "NTHRDS" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        if "ROOTPE" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        return value
#
# We need a set value until we full transition from perl
#

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        if "NTASKS" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        if "NTHRDS" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        if "ROOTPE" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        value = EnvBase.set_value(self, vid, value, subgroup, ignore_type)

