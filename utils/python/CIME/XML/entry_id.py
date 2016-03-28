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
            # loop through all the keys in valnode attributes
            for key,value in valnode.attrib.iteritems():
                # determine if key is in attributes dictionary
                if key in attributes:
                    if re.search(value, attributes[key]):
                        # set node value to valnode.text
                        node.set("value", valnode.text)
                        logger.debug("id %s value %s" % (node.get("id"), valnode.text))
                        return valnode.text

        # Fall back to default value
        value = self.get_optional_node("default_value", root=node)
        if value is not None and value.text is not None:
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

    def _set_value(self, node, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        if ignore_type:
            expect(type(value) is str, "Value must be type string if ignore_type is true")
            node.set("value",value)
        else:
            type_str = self._get_type_info(node)
            node.set("value", convert_to_string(value, type_str, vid))
        return value

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            return self._set_value(node, vid, value, subgroup, ignore_type)
        else:
            # Add item to GenericXML object in the lookup dictionary
            GenericXML.set_value(self,vid, value)  #***THIS IS NEW ***
            return None

    def get_value(self, vid, attribute={}, resolved=True, subgroup=None):
        """
        Get a value for entry with id attribute vid.
        or from the values field if the attribute argument is provided
        and matches
        """
        val = None
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            val = GenericXML.get_value(self, vid, resolved)
            return val

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

        # Return value as right type if we were able to fully resolve
        # otherwise, we have to leave as string.
        if "$" in val:
            return val
        else:
            type_str = self._get_type_info(node)
            return convert_to_type(val, type_str, vid)

    def get_values(self, vid, att, resolved=True):
        """
        If an entry includes a list of values return a dict matching each
        attribute to its associated value
        """
        values = {}
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return
        type_str = self._get_type_info(node)
        logger.debug("vid %s type %s"%(vid,type_str))

        valnodes = self.get_nodes("value", root=node)
        for valnode in valnodes:
            vatt = valnode.get(att)
            if resolved:
                values[vatt] = self.get_resolved_value(valnode.text)
            else:
                values[vatt] = valnode.text
            values[vatt] = convert_to_type(values[vatt], type_str, vid)

        return values

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
                elements.append(deepcopy(node))

        return elements

    def add_elements_by_group(self, srcobj, attlist, infile):
        """
        Add elements from srcdoc to self under the appropriate
        group element, entries to be added must have a child element
        <file> with value "infile"
        """
        # First get the list of entries in srcdoc with matching file children
        nodelist = srcobj.get_elements_from_child_content('file', infile)

        # For matchs found: Remove {<group>, <file>, <values>}
        # children from each entry and set the default value for the
        # new entries in self - putting the entries as children of
        # group elements in file $file
        for src_node in nodelist:
	    node  = deepcopy(src_node)
            gnode = src_node.find(".//group")
            gname = gnode.text

	    # If group with id=$gname does not exist in self.groups
	    # then create the group node and add it to infile file
            if gname not in self.groups.keys():
                newgroup = ET.Element("group")
                newgroup.set("id",gname)
                # initialize an empty list
                self.groups[gname] = newgroup
                self.add_child(newgroup)

	    # Set the default value, it may be determined by a regular
            # expression match to a dictionary value in attlist matching a
            # value attribute in node
            self.set_default_value(node, attlist)

            # Remove {<group>, <file>, <values>} from the entry element
            node = self.cleanupnode(node)

	    # Add the entry element to the group
            self.groups[gname].append(node)
            logger.debug ("Adding to group " + gname)

        return nodelist

    def cleanupnode(self, node):
        """
        Remove the <group>, <file>, <values> and <value> childnodes from node
        """
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
