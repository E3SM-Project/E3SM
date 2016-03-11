"""
Interface to the env_batch.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvBatch(EnvBase):

    def __init__(self, case_root=os.getcwd(), infile="env_batch.xml"):
        """
        initialize an object interface to file env_batch.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)

    def set_value(self, item, value, subgroup=None):
        val = None
        if subgroup is None:
            nodes = self.get_nodes("entry", {"id":item})
            for node in nodes:
                node.attrib["value"] = value
                val = value
        else:
            nodes = self.get_nodes("job",{"name":subgroup})
            for node in nodes:
                vnode = self.get_optional_node("entry", {"id":item}, root=node)
                if vnode is not None:
                    val = EnvBase.set_value(self, vnode, value)

        return val

    def get_value(self, item, attribute={}, resolved=True, subgroup=None):
        value = None
        nodes = self.get_nodes("entry",{"id":item})
        if len(nodes) > 0:
            value = {} # Not consistent with return values for other classes' get_value methods
            if subgroup is None:
                nodes = self.get_nodes("job")
            else:
                nodes = self.get_nodes("job", {"name":subgroup})

            for node in nodes:
                val = EnvBase.get_value(self, node,
                                        item, attribute, resolved)
                if val is not None:
                    value[node.attrib["name"]] = val

        return value

    def get_type_info(self, vid):
        nodes = self.get_nodes("entry",{"id":vid})
        type_info = None
        for node in nodes:
            new_type_info = self._get_type_info(node)
            if type_info is None:
                type_info = new_type_info
            else:
                expect( type_info == new_type_info,
                        "Inconsistent type_info for entry id=%s %s %s" % (vid, new_type_info, type_info))
        return type_info
