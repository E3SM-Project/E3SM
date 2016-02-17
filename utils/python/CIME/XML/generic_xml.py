"""
Common interface to XML files, this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""
from standard_module_setup import *
from xml.dom import minidom
from CIME.utils import expect, get_cime_root

class GenericXML(object):

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
            root = ET.Element("xml")
            root.set('version','1.0')
            self.root = ET.SubElement(root,"file")
            self.root.set('id',infile)
            self.tree = ET.ElementTree(root)

    def read(self, infile):
        """
        Read and parse an xml file into the object
        """
        logging.info("read: "+infile)
        self.tree = ET.parse(infile)
        self.root = self.tree.getroot()
        if ("version" in self.root.attrib):
            self.version = self.root.attrib["version"]
        else:
            self.version = "1.0"
        logging.info("File version is "+self.version)


    def write(self, outfile=None):
        """
        Write an xml file from data in self
        """
        if(outfile is None):
            outfile = self.filename
        logging.info("write: "+ outfile)
        xmlstr = ET.tostring(self.root)
        doc = minidom.parseString(xmlstr)
        with open(outfile,'w') as xmlout:
            doc.writexml(xmlout,addindent='  ')


    def get_node(self, nodename, attributes=None, root=None):
        """
        Get an xml element matching nodename with optional attributes
        """
        if(root is None):
            root = self.root
        nodes = []
        xpath = ".//"+nodename
        if(attributes is not None):
            keys = list(attributes.keys())
            # xml.etree has limited support for xpath and does not allow more than
            # one attribute in an xpath query so we query seperately for each attribute
            # and create a result with the intersection of those lists
            for key in keys:
                xpath = ".//%s[@%s=\'%s\']" % (nodename,key,attributes[key])
                newnodes = root.findall(xpath)
                if(not nodes):
                    nodes = newnodes
                else:
                    for node in nodes[:]:
                        if (node not in newnodes):
                            nodes.remove(node)
                if(not nodes):
                    return []
        else:
            nodes = root.findall(xpath)

        return nodes

    def add_child(self, node, root=None):
        """
        Add element node to self at root
        """
        if(root is None):
            root = self.root
        self.root.append(node)

    def get_value(self, item,resolved=True):
        """
        get_value is expected to be defined by the derived classes, if you get here it is an error.
        """
        logging.debug("Get Value for "+item)
        result = None
        if item in self.lookups.keys():
            result = self.lookups[item]
        if item in os.environ:
            result = os.environ.get(item)

        if (result is None):
            logging.info("No value available for item '%s'" % item)
        elif(resolved):
            result = self.get_resolved_value(result)

        return result

    def set_value(self,vid, value):
        valnodes = self.get_node(vid)
        if(valnodes):
            for node in valnodes:
                node.text = value
        else:
            self.lookups[vid] = value

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
        logging.debug("raw_value %s" % raw_value)
        reference_re = re.compile(r'\$(\w+)')
        env_ref_re   = re.compile(r'\$ENV\{(\w+)\}')
        item_data = raw_value

        if (item_data is None):
            return None

        for m in env_ref_re.finditer(item_data):
            logging.debug("look for "+item_data+ " in env")
            env_var = m.groups()[0]
            expect(env_var in os.environ, "Undefined env var '%s'" % env_var)
            item_data = item_data.replace(m.group(), os.environ[env_var])

        for m in reference_re.finditer(item_data):
            var = m.groups()[0]
            logging.debug("find: "+var)
            ref = self.get_value(var)
            if(ref is not None):
                logging.debug("resolve: "+ref)
                item_data = item_data.replace(m.group(), self.get_resolved_value(ref))

        return item_data
