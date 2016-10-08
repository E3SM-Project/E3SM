"""
Common interface to XML files which follow the entry id format,
this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""
from CIME.XML.standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type
from CIME.XML.generic_xml import GenericXML

from copy import deepcopy

logger = logging.getLogger(__name__)

class EntryID(GenericXML):

    def __init__(self, infile=None):
        GenericXML.__init__(self, infile)
        self.groups={}

    def get_default_value(self, node, attributes=None):
        """
        Set the value of an entry to the default value for that entry
        """
        value = self._get_value_match(node, attributes)
        if value is None:
            # Fall back to default value
            val_node = self.get_optional_node("default_value", root=node)
            if val_node is not None:
                value = val_node.text
        else:
            logger.debug("node is %s value is %s" % (node.get("id"), value))

        if value is None:
            logger.debug("For vid %s value is none"%node.get("id"))
            value = ""

        return value

    def set_default_value(self, vid, val):
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            default_node = self.get_optional_node("default_value", root=node)
            if default_node is not None:
                default_node.text = val
            else:
                logger.warn("Called set_default_value on a node without default_value field")

        return val

    def get_value_match(self, vid, attributes=None):
        # Handle this case:
        # <entry id ...>
        #  <values>
        #   <value A="a1">X</value>
        #   <value A="a2">Y</value>
        #   <value A="a3" B="b1">Z</value>
        #  </values>
        # </entry>

        node = self.get_optional_node("entry", {"id":vid})
        value = None
        if node is not None:
            value = self._get_value_match(node, attributes)
        return value

    def _get_value_match(self, node, attributes=None):
        '''
        Note that the component class has a specific version of this function
        '''
        match_value = None
        match_max = -1
        match_count = 0
        expect(node is not None," Empty node in _get_value_match")

        values = self.get_optional_node("values", root=node)
        if values is None:
            return

        for valnode in self.get_nodes("value", root=values):
            # loop through all the keys in valnode (value nodes) attributes
            if len(valnode.attrib) == 0:
                match_count = 0
            else:
                for key,value in valnode.attrib.iteritems():
                    # determine if key is in attributes dictionary
                    match_count = 0
                    if attributes is not None and key in attributes:
                        if re.search(value, attributes[key]):
                            logger.debug("Value %s and key %s match with value %s"%(value, key, attributes[key]))
                            match_count += 1
                        else:
                            match_count = -1
                            break

            if match_count > match_max:
                match_max = match_count
                match_value = valnode.text
            elif match_count == match_max:
                logger.debug("Ambiguous match for node '%s' for attributes '%s', falling back to order precedence" %
                             (node.attrib["id"], attributes))

        return match_value

    def get_node_element_info(self, vid, element_name):
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return None
        else:
            return self._get_node_element_info(node, element_name)

    def _get_node_element_info(self, node, element_name):
        element_node = self.get_optional_node(element_name, root=node)
        if element_node is not None:
            return element_node.text
        return None

    def _get_type_info(self, node):
        val = self._get_node_element_info(node, "type")
        if val is None:
            return "char"
        return val

    def get_type_info(self, vid):
        val = self.get_node_element_info(vid, "type")
        if val is None:
            return "char"
        return val

    def _get_default(self, node):
        return self._get_node_element_info(node, "default_value")

    # Get description , expect child with tag "description" for parent node
    def _get_description (self, node):
        return self._get_node_element_info(node, "desc")

    # Get group , expect node with tag "group"
    # entry id nodes are children of group nodes
    def _get_group (self, node):
        return self._get_node_element_info(node, "group")

    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        expect(subgroup is None, "Subgroup not supported")
        valid_values = self.get_optional_node("valid_values", root=node)
        if ignore_type:
            expect(type(value) is str, "Value must be type string if ignore_type is true")
            str_value = value
        else:
            type_str = self._get_type_info(node)
            str_value = convert_to_string(value, type_str, vid)
        if valid_values is not None and valid_values.text is not None and not str_value.startswith('$'):
            vvlist = [item.lstrip() for item in valid_values.text.split(',')]
            expect(str_value in vvlist, "Did not find %s in valid values:%s"%(value, vvlist))
        node.set("value", str_value)

        return value

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        val = None
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            val = self._set_value(node, value, vid, subgroup, ignore_type)
        return val

    def get_values(self, vid, attribute=None, resolved=True, subgroup=None):
        """
        Same functionality as get_value but it returns a list, if the
        value in xml contains commas the list have multiple elements split on
        commas
        """
        results = []
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return results
        str_result = self._get_value(node, attribute=attribute, resolved=resolved, subgroup=subgroup)
        str_results = str_result.split(',')
        for result in str_results:
            # Return value as right type if we were able to fully resolve
            # otherwise, we have to leave as string.
            if "$" in result:
                results.append(result)
            else:
                type_str = self._get_type_info(node)
                results.append( convert_to_type(result, type_str, vid))
        return results

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None):
        """
        Get a value for entry with id attribute vid.
        or from the values field if the attribute argument is provided
        and matches
        """
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return
        val = self._get_value(node, attribute=attribute, resolved=resolved, subgroup=subgroup)
        # Return value as right type if we were able to fully resolve
        # otherwise, we have to leave as string.
        if "$" in val:
            return val
        else:
            type_str = self._get_type_info(node)
            return convert_to_type(val, type_str, vid)

    def _get_value(self, node, attribute=None, resolved=True, subgroup=None):
        """
        internal get_value, does not convert to type
        """
        logger.debug("(_get_value) (%s, %s, %s)" % (attribute, resolved, subgroup))
        val = None
        if node is None:
            logger.debug("No node")
            return val

        logger.debug("Found node %s with attributes %s" , node.tag , node.attrib)
        if attribute:
            valnodes = self.get_optional_node("value", attribute, root=node)
            if valnodes is not None:
                val = valnodes.text
        elif node.get("value") is not None:
            val = node.get("value")
        else:
            val = self.get_default_value(node)

        if resolved:
            val = self.get_resolved_value(val)

        return val

    def get_full_records(self, item, attribute=None, resolved=True, subgroup=None): # (self, vid, att, resolved=True , subgroup=None ):
        """
        If an entry includes a list of values return a list of dict matching each
        attribute to its associated value and group
        """
        logger.debug("(get_full_records) Input values: %s , %s , %s , %s , %s" ,  self.__class__.__name__ , item, attribute, resolved, subgroup)

        nodes   = [] # List of identified xml elements
        results = [] # List of identified parameters

        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']

        groups = self.get_nodes("group")

        if (len(groups) == 0) :
            groups = [None]

        for group in groups :

            if item :
                nodes = self.get_nodes("entry",{"id" : item} , root=group)
            else :
                # Return all nodes
                logger.debug("(get_full_records) Retrieving all parameter")
                nodes = self.get_nodes("entry" , root=group)

            if (len(nodes) == 0) :
                logger.debug("(get_full_records) Found no nodes for %s" , item)
            else :
                logger.debug("(get_full_records) Building return structure for %s nodes" , len(nodes))

            for node in nodes :
                logger.debug("(get_full_records) Node tag=%s attribute=%s" , node.tag , node.attrib )

                g       = self._get_group(group)
                val     = node.attrib['value']
                attr    = node.attrib['id']
                t       = self._get_type_info(node)
                desc    = self._get_description(node)
                default = self._get_default(node)
                file_   = None
                try :
                    file_   = self.filename
                except AttributeError:
                    logger.debug("(get_full_records) Can't call filename on %s (%s)" , self , self.__class__.__name__ )
                #t   =  super(EnvBase , self).get_type( node )
                v = { 'group' : g , 'attribute' : attr , 'value' : val , 'type' : t , 'description' : desc , 'default' : default , 'file' : file_}

                results.append(v)

        logger.debug("(get_full_records) Returning %s items" , len(results) )
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
            childnode = self.get_optional_node(childname,root=node)
            expect(childnode is not None,"No childname %s for id %s"%(childname,node.get("id")))
            content = childnode.text
            if content == childcontent:
                elements.append(deepcopy(node))

        return elements

    def add_elements_by_group(self, srcobj, attributes=None, infile=None):
        """
        Add elements from srcobj to self under the appropriate
        group element, entries to be added must have a child element
        <file> with value "infile"
        """
        if infile is None:
            infile = os.path.basename(self.filename)
        # First get the list of entries in srcobj with matching file children

        nodelist = srcobj.get_elements_from_child_content('file', infile)

        # For matchs found: Remove {<group>, <file>, <values>}
        # children from each entry and set the default value for the
        # new entries in self - putting the entries as children of
        # group elements in file $file
        for src_node in nodelist:
            node  = deepcopy(src_node)
            gnode = src_node.find(".//group")
            gname = gnode.text
            if gname is None:
                gname = "group_not_set"
	    # If group with id=$gname does not exist in self.groups
	    # then create the group node and add it to infile file
            if gname not in self.groups.keys():
                newgroup = ET.Element("group")
                newgroup.set("id",gname)
                # initialize an empty list
                self.groups[gname] = newgroup
                self.add_child(newgroup)

            # Remove {<group>, <file>, <values>} from the entry element
            node = self.cleanupnode(node)

	    # Add the entry element to the group
            self.groups[gname].append(node)

	    # Set the default value, it may be determined by a regular
            # expression match to a dictionary value in attributes matching a
            # value attribute in node
            value = srcobj.get_default_value(src_node, attributes)

            self._set_value(node,value)
            logger.debug ("Adding to group " + gname)

        return nodelist

    def cleanupnode(self, node):
        """
        Remove the <group>, <file>, <values> and <value> childnodes from node
        """
        fnode = node.find(".//file")
        node.remove(fnode)
        gnode = node.find(".//group")
        node.remove(gnode)
        dnode = node.find(".//default_value")
        node.remove(dnode)
        vnode = node.find(".//values")
        if vnode is not None:
            node.remove(vnode)
        return node

    def compare_xml(self, other):
        xmldiffs = {}
        f1nodes = self.get_nodes("entry")
        for node in f1nodes:
            vid = node.get("id")
            f2val = other.get_value(vid, resolved=False)
            if f2val is not None:
                f1val = self.get_value(vid, resolved=False)
                if f2val != f1val:
                    xmldiffs[vid] = [f1val, f2val]

        return xmldiffs

    def __iter__(self):
        for node in self.get_nodes("entry"):
            vid = node.get("id")
            yield vid, self.get_value(vid)
