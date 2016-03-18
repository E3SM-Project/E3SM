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
        # allow the user to set all instances of item if subgroup is not provided
        if subgroup is None:
            nodes = self.get_nodes("entry", {"id":item})
            for node in nodes:
                self._set_value(node, item, value)
                val = value
        else:
            nodes = self.get_nodes("job",{"name":subgroup})
            for node in nodes:
                vnode = self.get_optional_node("entry", {"id":item}, root=node)
                if vnode is not None:
                    val = EnvBase.set_value(self, vnode, value)

        return val

    # def get_value(self, item, attribute={}, resolved=True, subgroup=None):
#         value = None
#         logger.debug("Value=%s additional Attribute:%s" , item , attribute)
#         nodes = self.get_nodes("entry",{"id":item})
#         logger.debug("Nodes:%s" , nodes)
#         if len(nodes) > 0:
#             value = {} # Not consistent with return values for other classes' get_value methods
#             if subgroup is None:
#                 nodes = self.get_nodes("job")
#             else:
#                 nodes = self.get_nodes("job", {"name":subgroup})
#
#             for node in nodes:
#                 logger.debug("Arguments for EnvBase.get_value: self=%s node=%s item=%s attribute=%s resolved=%s" , self, node , item , attribute , resolved)
#                 # Call parent method
#                 val = super(EnvBase , self).get_value( node,
#                                         item, attribute, resolved , node.attrib[""])
#                 # val = EnvBase.get_value(self, node,
# #                                         item, attribute, resolved)
#                 if val is not None:
#                     value[node.attrib["name"]] = val
#
#         return value

    def get_values(self, item, attribute={}, resolved=True, subgroup=None):
        """Returns the value as a string of the first xml element with item as attribute value. 
        <elememt_name attribute='attribute_value>value</element_name>"""
                
        nodes   = [] # List of identified xml elements  
        results = [] # List of identified parameters 
        logger.debug("Get node with attribute value: %s" , item)
       
        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        for job in self.get_nodes("job") :
            
            group = job.attrib['name']
            
            if item :
                nodes = self.get_nodes("*",{"id" : item} , root=job)
            else :
               # Return all nodes
               logger.debug("Retrieving all parameter")
               nodes = self.get_nodes("*",{"id" : None } , root=job)
          
           
            # Return value for first occurence of node with attribute value = item
            for node in nodes:
                t   =  super(EnvBase , self).get_type( node )
                val = { 'group' : group , 'attribute' : item , 'value' : node.attrib["value"] , 'type' : t , 'description' : self.get_description }
                logger.debug("Found node with value for %s = %s" , item , val)   
                results.append(val)
        
        return results

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
