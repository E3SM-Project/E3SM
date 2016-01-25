"""
Common interface to XML files which follow the entry id format, 
this is an abstract class and is expected to 
be used by other XML interface modules and not directly. 
"""

import xml.etree.ElementTree as ET
import os.path
import logging
import re
from GenericXML import GenericXML
from CIME.utils import expect

class EntryID(GenericXML):
    def __init__(self, infile=None):
        GenericXML.__init__(self,infile)
        
    def set_default_value(self, vid, attributes=None):
        """ Set the value of an entry to the default value for that entry
        vid can be an xml node pointer or a string identifier of a node """
        value = None
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.get_node("entry",{"id":vid})
            if(nodes is None):
                return
            expect(len(nodes) == 1,"More than one match found for id " + vid)
            node = nodes[0]
        valnodes = self.get_node("value",root=node)
        if(valnodes is not None):
            for valnode in valnodes:
                for att in valnode.attributes:
                    if(att.key in attributes):
                        if(re.search(attributes[att.key],att.text)):
                            value = valnode.text
                            logging.info("id %s value %s" % (vid, valnode.text))
        if(value is None):
            value = self.get_node("default_value",root=node)
        if(value is not None):
            node.set("value",value)
        return value[0].text

    def set_value(self, vid, value):
        val = None
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.get_node("entry",{"id":vid})
            if(nodes is None):
                return
            expect(len(nodes) == 0,"More than one match found for id " + vid)
            node = nodes[0] 
        if(node is not None):
            val = value
            node.set("value",value)
        return val

    def get_value(self, vid, attribute=None):
        val = None
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.get_node("entry",{"id":vid})
            if(nodes is None):
                return
            node = nodes[0]

        if(attribute is not None):
            valnodes = self.get_node("value",attribute)
            if(valnodes is not None and len(valnodes) == 1):
                val = valnodes[0].text
        elif(node.get("value") is not None):
            val = node.get("value")
        else:
            val = self.set_default_value(vid)
        return val

    def get_values(self,vid,att):
        values = {}
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.get_node("entry",{"id":vid})
            if(nodes is None):
                return
            node = nodes[0]
        valnodes = self.get_node("value",node=node)
        if(valnodes is not None):
            for valnode in valnodes:
                vatt = valnode.attrib(att)
                values[vatt] = valnode.text
        return values

    
