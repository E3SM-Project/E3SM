"""
Common interface to XML files, this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""

import xml.etree.ElementTree as ET
import sys, os, logging, re, doctest

_LIB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.append(_LIB_DIR)

import cime_util
from cime_util import expect

###############################################################################
class GenericXML:
###############################################################################

    ###########################################################################
    def __init__(self, infile=None):
    ###########################################################################
        """
        Initialize an object
        """
        self.tree = None
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
            self.root = ET.ElementTree(infile)
            self.root.ET.set('version','1.0')

    ###########################################################################
    def read(self,infile):
    ###########################################################################
        """
        Read and parse an xml file into the object
        """
        logging.debug("read: "+infile)
        self.tree = ET.parse(infile)
        self.root = self.tree.getroot()

    ###########################################################################
    def write(self, cimeroot, infile):
    ###########################################################################
        """
        Write an xml file from data in self
        """
        logging.debug("write: "+ infile)
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

    ###########################################################################
    def get_node(self, nodename, attributes=None):
    ###########################################################################
        """
        Get an xml element matching nodename with optional attributes
        """
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
            return nodes[0]
        elif(findcnt > 1):
            return nodes
        else:
            return None

    ###########################################################################
    def add_child(self, node):
    ###########################################################################
        """
        Add element node to self at root
        """
        self.tree.insert(self.root, node)

    ###########################################################################
    def get_value(self, item):
    ###########################################################################
        """
        TODO
        """
        expect(False, "Not implemented")

    ###########################################################################
    def get_resolved_value(self, raw_value):
    ###########################################################################
        """
        A value in the xml file may contain references to other xml
        variables or to environment variables. These are refered to in
        the perl style with $name and $ENV{name}.

        >>> obj = GenericXML()
        >>> os.environ["FOO"] = "BAR"
        >>> os.environ["BAZ"] = "BARF"
        >>> obj.get_resolved_value("one $ENV{FOO} two $ENV{BAZ} three")
        'one BAR two BARF three'
        """
        reference_re = re.compile(r'\$(\w+)')
        env_ref_re   = re.compile(r'\$ENV\{(\w+)\}')

        item_data = raw_value
        for m in env_ref_re.finditer(item_data):
            env_var = m.groups()[0]
            expect(env_var in os.environ, "Undefined env var '%s'" % env_var)
            item_data = item_data.replace(m.group(), os.environ[env_var])

        for m in reference_re.finditer(item_data):
            ref = m.groups()[0]
            item_data = item_data.replace(m.group(), self.get_resolved_value(ref))

        return item_data
