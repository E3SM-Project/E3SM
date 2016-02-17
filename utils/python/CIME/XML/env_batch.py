"""
Interface to the env_batch.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

class EnvBatch(EnvBase):
    def __init__(self, case_root=os.getcwd(), infile="env_batch.xml"):
        """
        initialize an object interface to file env_batch.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)

    def set_value(self, item, value, subgroup=None):
        val = None
        if(subgroup is None):
            nodes = self.get_node("entry",{"id":item})
            for node in nodes:
                val = EnvBase.set_value(self, node, value)
        else:
            nodes = self.get_node("job",{"name":subgroup})
            for node in nodes:
                vnode = self.get_node("entry",{"id":item},root=node)
                if( len(vnode)>0 ):
                    val = EnvBase.set_value(self,vnode[0],value)

        return val

    def get_value(self, item, attribute={}, resolved=True, subgroup=None):
        value = {}
        nodes = self.get_node("entry",{"id":item})
        if(len(nodes)>0):
            if(subgroup is None):
                nodes = self.get_node("job")
            else:
                nodes = self.get_node("job",{"name":subgroup})
            for node in nodes:
                val = EnvBase.get_value(self, node,
                                        item,attribute,resolved)
                if(val):
                    value[node.attrib["name"]] = val
            return value







