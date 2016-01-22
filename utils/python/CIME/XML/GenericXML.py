"""
Common interface to XML files, this is an abstract class and is expected to 
be used by other XML interface modules and not directly. 
"""

import xml.etree.ElementTree as ET
import os.path
import logging

class GenericXML:
    def __init__(self, infile=None):
        """ Initialize an object """
        if(infile == None):
            # if file is not defined just return
            self.filename = None
            return 
        if(os.path.isfile(infile) and os.access(infile, os.R_OK)):
            # If file is defined and exists, read it
            self.filename = infile
            self.read(infile)
        else:
            # if file does not exist create a root xml element 
            # and set it's id to file
            self.filename = infile
            self.root = ET.Element(infile)
            self.root.set('version','1.0')
            
    def read(self,infile):
        """ Read and parse an xml file into the object """
        logging.info("read: "+infile)
        self.tree = ET.parse(infile)
        self.root = self.tree.getroot()

    def write(self, cimeroot, infile):
        """ Write an xml file from data in self """
        logging.info("write: "+ infile)
#        styledoc = os.path.join(cimeroot,"cime_config","xml_schemas","case_xml.xsl")
#        stylesheet = libxslt.parseStylesheetDoc(styledoc)
#        result = stylesheet.applyStyleSheet(self.root, None)
#        style.saveResultToFilename(file, result, 0)
#        style.freeStylesheet()
#        result.freeDoc()
        if(infile != None):
            self.tree.write(infile)
        else:
            self.tree.write(self.filename)
    
    def GetNode(self, nodename, attributes=None):
        """ Get an xml element matching nodename with optional attributes """
        xpath = ".//"+nodename
        if(attributes):
            keys = list(attributes.keys())
            cnt = 0
            for key in keys:
                if(cnt == 0):
                    xpath += "["
                else:
                    xpath += " and "
                xpath += "@%s=\'%s\'" % (key,attributes[key])
                cnt=cnt+1
            xpath += "]"
        nodes = self.root.findall(xpath)
        findcnt = len(nodes)
        if(findcnt == 1):
            return node[0]
        elif(findcnt > 1):
            return nodes
        else:
            return None

    def addChild(self, node):
        """ Add element node to self at root """
        self.tree.insert(self.root, node)



