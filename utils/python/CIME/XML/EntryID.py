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

class EntryID(GenericXML):
    def __init__(self, infile=None):
        GenericXML.__init__(self,infile)
        
    def SetDefaultValue(self, vid, attributes):
        """ Set the value of an entry to the default value for that entry
        vid can be an xml node pointer or a string identifier of a node """
        value = None
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.GetNode("entry",{"id":vid})
            if(nodes is None):
                return
            expect(len(nodes) == 0,"More than one match found for id " + vid)
            node = nodes[0]
        valnodes = self.GetNode("value",root=node)
        if(valnodes is not None):
            for valnode in valnodes:
                for att in valnode.attributes:
                    if(att.key in attributes):
                        if(re.search(attributes[att.key],att.text)):
                            value = valnode.text
                            logging.info("id %s value %s" % (vid, valnode.text))
        if(value is None):
            value = self.GetNode("default_value",root=node)
        if(value is not None):
            node.set("value",value)
        return value

    def SetValue(self, vid, value):
        val = None
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.GetNode("entry",{"id":vid})
            if(nodes is None):
                return
            expect(len(nodes) == 0,"More than one match found for id " + vid)
            node = nodes[0] 
        if(node is not None):
            val = value
            node.set("value",value)
        return val

    def GetValue(self, vid, attribute=None):
        val = None
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.GetNode("entry",{"id":vid})
            if(nodes is None):
                return
            node = nodes[0]
        if(attribute is not None):
            valnodes = self.GetNode("value",attribute)
            if(valnodes is not None and len(valnodes) == 1):
                val = valnodes[0].text
        elif(node.attrib("value")):
            val = node.attrib("value")
        else:
            self.SetDefaultValue(vid,value)

    def GetValues(self,vid,att):
        values = {}
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self.GetNode("entry",{"id":vid})
            if(nodes is None):
                return
            node = nodes[0]
        valnodes = self.GetNode("value",node=node)
        if(valnodes is not None):
            for valnode in valnodes:
                vatt = valnode.attrib(att)
                values[vatt] = valnode.text
        return values

    
