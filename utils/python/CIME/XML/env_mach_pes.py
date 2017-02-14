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

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None, pes_per_node=None): # pylint: disable=arguments-differ
        value = EnvBase.get_value(self, vid, attribute, resolved, subgroup)
        if "NTASKS" in vid or "ROOTPE" in vid and pes_per_node is None:
            pes_per_node = self.get_value("PES_PER_NODE")

            if "NTASKS" in vid and value < 0:
                value = -1*value*pes_per_node
            if "ROOTPE" in vid and value < 0:
                value = -1*value*pes_per_node
        return value

    def set_value(self, vid, value, subgroup=None, ignore_type=False, pes_per_node=None): # pylint: disable=arguments-differ
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        val = None
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            val = self._set_value(node, value, vid, subgroup, ignore_type, pes_per_node=pes_per_node)
        return val



    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False, pes_per_node=None): # pylint: disable=arguments-differ
        if vid is None:
            vid = node.get("id")

        if "NTASKS" in vid or "ROOTPE" in vid and pes_per_node is None:
            pes_per_node = self.get_value("PES_PER_NODE")

        if "NTASKS" in vid and value < 0:
            value = -1*value*pes_per_node
        if "ROOTPE" in vid and value < 0:
            value = -1*value*pes_per_node
        val = EnvBase._set_value(self, node, value, vid, subgroup, ignore_type)
        return val
