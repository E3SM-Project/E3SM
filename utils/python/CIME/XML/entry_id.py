"""
Common interface to XML files which follow the entry id format,
this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""
from standard_module_setup import *
from CIME.utils import expect
from generic_xml import GenericXML

class EntryID(GenericXML):

    def __init__(self, infile=None):
        GenericXML.__init__(self, infile)
        self.groups={}

    def set_default_value(self, vid, attributes=None):
        """
        Set the value of an entry to the default value for that entry
        vid can be an xml node pointer or a string identifier of a node
        """
        value = None
        if (type(vid) != type(str())):
            node = vid
            vid = node.attrib["id"]
        else:
            nodes = self.get_node("entry", {"id":vid})
            if (nodes is None):
                return
            expect(len(nodes) == 1, "More than one match found for id " + vid)
            node = nodes[0]

        valnodes = self.get_node("value",root=node)
        if (valnodes is not None):
            for valnode in valnodes:
                for att in valnode.attributes:
                    if (att.key in attributes):
                        if (re.search(attributes[att.key],att.text)):
                            value = valnode.text
                            logging.info("id %s value %s" % (vid, valnode.text))

        if (value is None):
            value = self.get_node("default_value", root=node)
        if (value is not None):
            node.set("value", value[0].text)
            return value[0].text

    def set_value(self, vid, value,subgroup=None):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        val = None
        if (type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.get_node("entry", {"id":vid})
            if (not nodes):
                return
            expect(len(nodes) == 1, "More than one match found for id " + vid)
            node = nodes[0]

        if (node is not None):
            val = value
            node.set("value", value)

        return val

    def get_value(self, vid, attribute={}, resolved=True,subgroup=None):
        """
        get a value for entry with id attribute vid.
        or from the values field if the attribute argument is provided
        and matches
        """
        val = None
        if (type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.get_node("entry", {"id":vid})
            if (len(nodes) == 0):
                val = GenericXML.get_value(self, vid, resolved)
                return val
            else:
                node = nodes[0]

        if (attribute):
            valnodes = self.get_node("value", attribute,root=node)
            if (valnodes is not None and len(valnodes) == 1):
                val = valnodes[0].text
        elif (node.get("value") is not None):
            val = node.get("value")
        else:
            val = self.set_default_value(vid)

        if (val is None):
            # if all else fails
            val = GenericXML.get_value(self,vid,resolved)

        if (resolved):
            val = self.get_resolved_value(val)

        return val

    def get_values(self, vid, att, resolved=True):
        """
        If an entry includes a list of values return a dict matching each
        attribute to its associated value
        """
        values = {}
        if (type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.get_node("entry", {"id":vid})
            if (nodes is None):
                return
            node = nodes[0]

        valnodes = self.get_node("value", root=node)
        if (valnodes is not None):
            for valnode in valnodes:
                vatt = valnode.attrib[att]
                if(resolved):
                    values[vatt] = self.get_resolved_value(valnode.text)
                else:
                    values[vatt] = valnode.text


        return values

    def get_child_content(self,id,childname):
        val = None
        node = self.get_node("entry",{"id":id})
        if(node):
            childmatch  = self.get_node(childname, root=node[0])
            if(len(childmatch) == 1):
                val = childmatch[0].text
        return val

    def get_elements_from_child_content(self, childname, childcontent):
        nodes = self.get_node("entry")
        elements = []
        for node in nodes:
            childnodes = self.get_node(childname,root=node)
            expect(len(childnodes)==1,"Unexpected number of matchs for %s in %s"%(childname, node.get("id")))
            content = childnodes[0].text
            if(content == childcontent):
                elements.append( node)

        return elements

    def add_elements_by_group(self, srcobj, attlist, infile):
        nodelist = srcobj.get_elements_from_child_content('file',infile)
        for node in nodelist:
            gnode = node.find(".//group")
            gname = gnode.text
            if(gname not in self.groups.keys()):
                newgroup = ET.Element("group")
                newgroup.set("id",gname)
                self.add_child(newgroup)
                self.groups[gname] = newgroup
            self.set_default_value(node, attlist)
            node = self.cleanupnode(node)
            self.groups[gname].append(node)
            logging.info ("Adding to group "+gname)
        return nodelist

    def cleanupnode(self,node):
        fnode = node.find(".//file")
        gnode = node.find(".//group")
        node.remove(fnode)
        node.remove(gnode)
        vnode = node.find(".//values")
        if(vnode):
            node.remove(vnode)
        return node
