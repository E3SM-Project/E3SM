"""
Common interface to XML files which follow the entry id format,
this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""
from standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type
from generic_xml import GenericXML

logger = logging.getLogger(__name__)

class EntryID(GenericXML):

    def __init__(self, infile=None):
        GenericXML.__init__(self, infile)
        self.groups={}

    def set_default_value(self, node, attributes=None):
        """
        Set the value of an entry to the default value for that entry
        """
        # Hangle this case:
        # <entry id ...>
        #  <values>
        #   <value A="a1">X</value>
        #   <value A="a2">Y</value>
        #   <value A="a3">Z</value>
        #  </values>
        # </entry>
        valnodes = self.get_nodes("value", root=node)
        for valnode in valnodes:
            for att in valnode.attributes:
                if att.key in attributes:
                    if re.search(attributes[att.key], att.text):
                        node.set("value", valnode.text)
                        logger.debug("id %s value %s" % (node.attrib["id"], valnode.text))
                        return valnode.text

        # Fall back to default value
        value = self.get_optional_node("default_value", root=node)
        if value is not None:
            node.set("value", value.text)
            return value.text

    def _get_type_info(self, node):
        type_node = self.get_optional_node("type", root=node)
        if type_node is not None:
            return type_node.text
        else:
            # Default to string
            return "char"

    def get_type_info(self, vid):
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return None
        else:
            return self._get_type_info(node)
    
    def get_default(self, node):
        default = self.get_optional_node("default_value", root=node)
        if default is not None:
            return default.text
        else:
            return None
        
    
    # Get type description , expect child with tag "type" for node         
    def get_type (self, node):
        type_node = self.get_optional_node("type", root=node)
        if type_node is not None:
            return type_node.text
        else:
            # Default to string
            return "char"
            
    # Get description , expect child with tag "description" for parent node
    def get_description (self, node):
        type_node = self.get_optional_node("desc", root=node)
        if type_node is not None:
            return type_node.text
        else:
            return None

    # Get group , expect child with tag "group" for parent node        
    def get_group (self, node):
        type_node = self.get_optional_node("group", root=node)
        if type_node is not None:
            return type_node.text
        else:
            # Default to None
            return None

    def _set_value(self, node, vid, value, subgroup=None):
        type_str = self._get_type_info(node)
        node.set("value", convert_to_string(value, type_str, vid))
        return value

    def set_value(self, vid, value, subgroup=None):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            return self._set_value(node, vid, value, subgroup)
        else:
            return None

    def get_value(self, vid, attribute={}, resolved=True, subgroup=None):
        """
        get a value for entry with id attribute vid.
        or from the values field if the attribute argument is provided
        and matches
        """
        
        logger.debug("Get Value")
        val = None
        node = self.get_optional_node("entry", {"id":vid})

        if node is None:
            logger.debug("No node")
            val = GenericXML.get_value(self, vid, resolved)
            return val

        logger.debug("Found node %s with attributes %s" , node.tag , node.attrib)
        if attribute:
            valnodes = self.get_optional_node("value", attribute, root=node)
            if valnodes is not None:
                val = valnodes.text
        elif node.get("value") is not None:
            val = node.get("value")
        else:
            val = self.set_default_value(node)

        if val is None:
            # if all else fails
            val = GenericXML.get_value(self,vid,resolved)

        if resolved:
            val = self.get_resolved_value(val)

        # Return value as right type
        type_str = self._get_type_info(node)
        return convert_to_type(val, type_str, vid)

    def get_values(self, item, attribute={}, resolved=True, subgroup=None): # (self, vid, att, resolved=True , subgroup=None ):
     
        """
        If an entry includes a list of values return a list of dict matching each
        attribute to its associated value and group
        """
        
        nodes   = [] # List of identified xml elements  
        results = [] # List of identified parameters 
        logger.debug("Here Get node with attribute value: %s" , item)
       
        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        
        if item :
            nodes = self.get_nodes("entry",{"id" : item})
        else :
           # Return all nodes
           logger.debug("Retrieving all parameter")
           nodes = self.get_nodes("entry")
        
        if (len(nodes) == 0) :
            logger.debug("Found no nodes for %s" , item)
        else :
            logger.debug("Building return structure for %s nodes" , len(nodes))
        
        for node in nodes :
            logger.debug("Node tag=%s attribute=%s" , node.tag , node.attrib )
            
            group   = self.get_group(node)
            val     = node.attrib['value']
            attr    = node.attrib['id']
            t       = self.get_type(node)
            desc    = self.get_description(node)
            default = self.get_default(node)
            file    = None
            try :     
                file    = self.filename
            except AttributeError:
                logger.debug("Can't call filename on %s (%s)" , self , self.__class__.__name__ )      
            #t   =  super(EnvBase , self).get_type( node )
            v = { 'group' : group , 'attribute' : attr , 'value' : val , 'type' : t , 'description' : desc , 'default' : default , 'file' : file}
            logger.debug("Found node with value for %s = %s" , item , v )   
            results.append(v)
        
        return results
       
   
 

    def get_child_content(self, vid, childname):
        val = None
        node = self.get_optional_node("entry", {"id" : vid})
        if node is not None:
            childmatch  = self.get_optional_node(childname, root=node[0])
            if childmatch is not None:
                val = childmatch.text

        return val

    def get_elements_from_child_content(self, childname, childcontent):
        nodes = self.get_nodes("entry")
        elements = []
        for node in nodes:
            childnode = self.get_node(childname,root=node)
            content = childnode.text
            if content == childcontent:
                elements.append(node)

        return elements

    def add_elements_by_group(self, srcobj, attlist, infile):
        nodelist = srcobj.get_elements_from_child_content('file', infile)
        for node in nodelist:
            gnode = node.find(".//group")
            gname = gnode.text
            if gname not in self.groups.keys():
                newgroup = ET.Element("group")
                newgroup.set("id",gname)
                self.add_child(newgroup)
                self.groups[gname] = newgroup

            self.set_default_value(node, attlist)
            node = self.cleanupnode(node)
            self.groups[gname].append(node)
            logger.debug ("Adding to group " + gname)

        return nodelist

    def cleanupnode(self, node):
        fnode = node.find(".//file")
        gnode = node.find(".//group")
        node.remove(fnode)
        node.remove(gnode)
        vnode = node.find(".//values")
        if vnode is not None:
            node.remove(vnode)
        return node

    def compare_xml(self, other):
        xmldiffs = {}
        f1nodes = self.get_nodes("entry")
        for node in f1nodes:
            vid = node.attrib["id"]
            f2val = other.get_value(vid, resolved=False)
            if f2val is not None:
                f1val = self.get_value(vid, resolved=False)
                if f2val != f1val:
                    xmldiffs[vid] = [f1val, f2val]

        return xmldiffs

