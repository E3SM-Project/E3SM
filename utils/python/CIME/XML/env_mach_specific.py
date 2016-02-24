"""
Interface to the env_mach_specific.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

class EnvMachSpecific(EnvBase):

    def __init__(self, caseroot=os.getcwd(), infile="env_mach_specific.xml"):
        """
        initialize an object interface to file env_mach_specific.xml in the case directory
        """
        EnvBase.__init__(self, caseroot, infile)

    
    def get_values(self, item="*", attribute={}, resolved=True, subgroup=None):
            
        value = {}  
        root  = self.tree.getroot()
        xpath = [ ".//" + item + "[@name]" ]
        xpathes = []
        
        # create search path for all attributes      
        if bool(attribute) :
            for key, attribute_value in attribute.iteritems():
                if attribute_value :
                    xpathes.append( ".//" + item + "[@" + key + "='" + attribute_value +"']" )
                else:
                    xpathes.append( ".//" + item + "[@" + key + "]" )
        else:
            xpathes.append(xpath)
        
        # Get all values and return results as dicitonary     
        for xpath in xpathes:   
            logging.debug("XPATH (env_mach_specific.py): "+ xpath)
            
            for node in root.findall(xpath):
                if node.attrib:
                    for a in node.attrib:
                        value[node.attrib[a]] = node.text
                else:
                  value[node.tag] = node.text  
                  
        if value:
            return value
           
        return None
       

        
        
    def get_value(self, item, attribute={}, resolved=True, subgroup=None):
        """Returns the value as a string of the first xml element with item as attribute value. 
        <elememt_name attribute='attribute_value>value</element_name>"""
                
        value = {}
        nodes = []
        
        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        if item :
            nodes = self.get_node("*",{"name" : item})
            
            
        for node in nodes:
            val = node.text    
            if val:
                return val
        
        return None
            