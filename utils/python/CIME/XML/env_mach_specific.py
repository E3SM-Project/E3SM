"""
Interface to the env_mach_specific.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.generic_xml import GenericXML
from CIME.XML.env_base import EnvBase

logger = logging.getLogger(__name__)

# Is not of type EntryID but can use functions from EntryID (e.g
# get_type) otherwise need to implement own functions and make GenericXML parent class
class EnvMachSpecific(EnvBase):

    def __init__(self, caseroot, infile="env_mach_specific.xml"):
        """
        initialize an object interface to file env_mach_specific.xml in the case directory
        """
        fullpath = infile if os.path.isabs(infile) else os.path.join(caseroot, infile)
        EnvBase.__init__(self, fullpath)

    def get_values(self, item, attribute={}, resolved=True, subgroup=None):
        """Returns the value as a string of the first xml element with item as attribute value.
        <elememt_name attribute='attribute_value>value</element_name>"""

        logger.debug("(get_values) Input values: %s , %s , %s , %s , %s" , self.__class__.__name__ , item, attribute, resolved, subgroup)

        nodes   = [] # List of identified xml elements
        results = [] # List of identified parameters

        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        if item :
            nodes = self.get_nodes("*",{"name" : item})
        else :
            # Return all nodes
            logger.debug("Retrieving all parameter")
            nodes = self.get_nodes("env")

        # Return value for first occurence of node with attribute value = item
        for node in nodes:

            group   = super(EnvMachSpecific, self)._get_group(node)
            val     = node.text
            attr    = node.attrib['name']
            t       = self._get_type(node)
            desc    = self._get_description(node)
            #default = super(EnvBase , self).get_default(node)
            default = self._get_default(node)
            filename    = self.filename

            #t   =  super(EnvBase , self).get_type( node )
            v = { 'group' : group , 'attribute' : attr , 'value' : val , 'type' : t , 'description' : desc , 'default' : default , 'file' : filename}
            logger.debug("Found node with value for %s = %s" , item , v )
            results.append(v)

        return results
