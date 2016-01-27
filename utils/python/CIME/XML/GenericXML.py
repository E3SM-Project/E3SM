"""
Common interface to XML files, this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""

import xml.etree.ElementTree as ET
import sys, os, logging, re, doctest

from CIME.utils import expect, get_cime_root

class GenericXML:

    def __init__(self, infile=None):
        """
        Initialize an object
        """
        self.tree = None
        self.lookups = {}
        self.lookups['CIMEROOT'] = get_cime_root()
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
        """
        Read and parse an xml file into the object
        """
        logging.debug("read: "+infile)
        self.tree = ET.parse(infile)
        self.root = self.tree.getroot()

    def write(self, cimeroot, infile):
        """
        Write an xml file from data in self
        """
        logging.debug("write: "+ infile)
        if(infile != None):
            self.tree.write(infile)
        else:
            self.tree.write(self.filename)
    
    def get_node(self, nodename, attributes=None, root=None):
        """
        Get an xml element matching nodename with optional attributes
        """
        if(root is None):
            root = self.root

        xpath = ".//"+nodename
        if(attributes is not None):
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
        logging.info("xpath = "+ xpath)
        nodes = root.findall(xpath)
        return nodes

    def add_child(self, node):
        """
        Add element node to self at root
        """
        self.tree.insert(self.root, node)

    def get_value(self, item):
        """
        get_value is expected to be defined by the derived classes, if you get here it is an error.
        """
        logging.debug("Get Value for "+item)
        if item in self.lookups.keys():
            return self.lookups[item] 
        if item in os.environ:
            return os.environ.get(item)
#        expect(False, "Not implemented")

    def get_resolved_value(self, raw_value):
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
        logging.debug("raw_value "+raw_value)
        reference_re = re.compile(r'\$(\w+)')
        env_ref_re   = re.compile(r'\$ENV\{(\w+)\}')
        item_data = raw_value
        logging.debug(" item_data "+item_data)
        for m in env_ref_re.finditer(item_data):
            logging.debug("look for "+item_data+ " in env")
            env_var = m.groups()[0]
            expect(env_var in os.environ, "Undefined env var '%s'" % env_var)
            item_data = item_data.replace(m.group(), os.environ[env_var])

        for m in reference_re.finditer(item_data):
            var = m.groups()[0]
            logging.debug("find: "+var)
            ref = self.get_value(var)
            logging.debug("resolve: "+ref)
            item_data = item_data.replace(m.group(), ref)

        return item_data
