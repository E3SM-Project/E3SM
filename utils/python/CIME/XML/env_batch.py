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
        if(subgroup is not None):
            node = self.get_node("job",{"name":subgroup})
            if(node):
                vnode = self.get_node(item,root=node[0])
                EnvBase.set_value(self, vnode, value)
        else:
            EnvBase.set_value(self,item,value)

    def get_value(self, item, attribute={}, resolved=True, subgroup=None):
        value = None
        if(subgroup is not None):
            node = self.get_node("job",{"name":subgroup})
            if(node):
                item = self.get_node(item,root=node[0])

        value = EnvBase.get_value(self, item, attribute, resolved, subgroup)
        return value
