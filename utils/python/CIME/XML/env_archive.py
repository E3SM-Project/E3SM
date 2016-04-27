"""
Interface to the env_archive.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from CIME.XML.generic_xml import GenericXML
from CIME.XML.archive import Archive
from CIME.XML.headers import Headers

logger = logging.getLogger(__name__)

class EnvArchive(GenericXML):

    def __init__(self, case_root=None, infile="env_archive.xml"):
        """
        initialize an object interface to file env_archive.xml in the case directory
        """
        logger.debug("Case_root = %s" , case_root)
        if case_root is None:
            case_root = os.getcwd()
        if os.path.isabs(infile):
            fullpath = infile
        else:
            fullpath = os.path.join(case_root, infile)
            logger.debug("Fullpath = %s" , fullpath)
        GenericXML.__init__(self, fullpath)
        if not os.path.isfile(infile):
            headerobj = Headers()
            headernode = headerobj.get_header_node(os.path.basename(fullpath))
            self.root.append(headernode)
            archive = Archive()
            self.root.append(archive.root)


    def get_values(self, item, attribute={}, resolved=True, subgroup=None):
        """Returns the value as a string of the first xml element with item as attribute value. 
        <elememt_name attribute='attribute_value>value</element_name>"""
                
        nodes   = [] # List of identified xml elements  
        results = [] # List of identified parameters 
        logger.debug("Get node with attribute value: %s" , item)
       
        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        if item :
            nodes = self.get_nodes("*",{"name" : item})
        else :
           # Return all nodes
           logger.debug("Retrieving all parameter")
           nodes = self.get_nodes("comp_archive_spec")
          
           
        # Return value for first occurence of node with attribute value = item
        for node in nodes:
            
            group   = node.attrib['name']
            rootdir = node.find('./rootdir').text
            
            for file in node.findall('./files/file_extension') :
                
            
                keep    = file.find('./keep_last_in_rundir')
            
                val     = keep.text
                attr    = rootdir + "/" + file.find('./subdir').text + "/" + file.get('regex_suffix')
                t       = self.get_type(node)
                desc    = self.get_description(keep)
                #default = super(EnvBase , self).get_default(node)
                default = self.get_default(keep)
                file    = self.filename
            
                #t   =  super(EnvBase , self).get_type( node )
                v = { 'group' : group , 'attribute' : attr , 'value' : val , 'type' : t , 'description' : desc , 'default' : default , 'file' : file}
                logger.debug("Found node with value for %s = %s" , item , v )   
                results.append(v)
            
            
            # val = { 'group' : None , 'attribute' : item , 'value' : node.text , 'type' : None , 'description' : None }
#             logger.debug("Found node with value for %s = %s" , item , val)
#             results.append(val)
        
        return results 
        
    def get_type(self, node=None) :
        return None    
        
    def get_description(self, node=None):
        
        description = None
        logger.warning("Node %s" , node.tag )
        if node.tag == "keep_last_in_rundir" :
            if node.get('description') :
                description = node.get('description')
            elif node.find('./description') :
                description = node.find('./description').text
            else:      
                description = "Keep last in rundir"
        
        return description
        
    def get_default(self, node=None):
        return None
        
    def get_archive():
        pass
        
    class Archive:
        """Class for a single archive component"""
        
        def __init__ (self):
            self.name = None # <comp_archive_spec name=" $NAME ">
            pass
            
        def name(self , name=None):
            return self.name
        
            