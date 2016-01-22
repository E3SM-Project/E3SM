"""
Common interface to XML files which follow the entry id format, 
this is an abstract class and is expected to 
be used by other XML interface modules and not directly. 
"""

import xml.etree.ElementTree as ET
import os.path
import logging

class EntryID(GenericXML):
    def __init__(self, infile=None):
        GenericXML.__init__(self,infile)
        
    def SetDefaultValue(self, vid, attributes):
        """ Set the value of an entry to the default value for that entry
        vid can be an xml node pointer or a string identifier of a node """
        if(type(vid) != type(str())):
            node = vid
            vid = node.attrib("id")
        else:
            nodes = self->GetNode("entry",{"id":vid})
            if(nodes is None):
                return
            expect(len(nodes) == 0,"More than one match found for id " + vid)
            node = nodes[0]
        valnode = self->GetNode("value",root=node);
