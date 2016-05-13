"""
Interface to the env_mach_pes.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvMachPes(EnvBase):

    def __init__(self, case_root=None, infile="env_mach_pes.xml"):
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

    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False):
        if vid is None:
            vid = node.get("id")

        if "NTASKS" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        if "NTHRDS" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        if "ROOTPE" in vid and value < 0:
            value = -1*value*self.get_value("PES_PER_NODE")
        val = EnvBase._set_value(self, node, value, vid, subgroup, ignore_type)
        return val
