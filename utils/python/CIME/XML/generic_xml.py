"""
Common interface to XML files, this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""
from standard_module_setup import *
from xml.dom import minidom
from CIME.utils import expect, get_cime_root

logger = logging.getLogger(__name__)

class GenericXML(object):

    def __init__(self, infile=None):
        """
        Initialize an object
        """
        self.tree = None

        # Hold arbitary values. In create_newcase we may set values
        # for xml files that haven't been created yet. We need a place
        # to store them until we are ready to create the file. At file
        # creation we get the values for those fields from this lookup
        # table and then remove the entry. This was what I came up
        # with in the perl anyway and I think that we still need it
        # here.
        self.lookups = {}
        self.lookups['CIMEROOT'] = get_cime_root()

        if infile == None:
            # if file is not defined just return
            self.filename = None
            return
        if os.path.isfile(infile) and os.access(infile, os.R_OK):
            # If file is defined and exists, read it
            self.filename = infile
            self.read(infile)
        else:
            # if file does not exist create a root xml element
            # and set it's id to file
            self.filename = infile
            root = ET.Element("xml")
            root.set("version", "1.0")
            self.root = ET.SubElement(root, "file")
            self.root.set("id", infile)
            self.tree = ET.ElementTree(root)

    def read(self, infile):
        """
        Read and parse an xml file into the object
        """
        logger.debug("read: " + infile)
        self.tree = ET.parse(infile)
        self.root = self.tree.getroot()
        self.version = self.root.get("version")
        self.version = "1.0" if self.version is None else self.version
        logger.debug("File version is "+self.version)

    def write(self, outfile=None):
        """
        Write an xml file from data in self
        """
        if outfile is None:
            outfile = self.filename

        logger.debug("write: " + outfile)
        xmlstr = ET.tostring(self.root)
        doc = minidom.parseString(xmlstr)
        with open(outfile,'w') as xmlout:
            doc.writexml(xmlout,addindent='  ')

    def get_node(self, nodename, attributes=None, root=None):
        """
        Get an xml element matching nodename with optional attributes.

        Error unless exactly one match.
        """
        nodes = self.get_nodes(nodename, attributes=attributes, root=root)
        expect(len(nodes) == 1, "Incorrect number of matches, %d, for nodename '%s' and attrs '%s' in file '%s'" %
               (len(nodes), nodename, attributes, self.filename))
        return nodes[0]

    def get_optional_node(self, nodename, attributes=None, root=None):
        """
        Get an xml element matching nodename with optional attributes.

        Return None if no match.
        """
        nodes = self.get_nodes(nodename, attributes=attributes, root=root)

        expect(len(nodes) <= 1, "Multiple matches for nodename '%s' and attrs '%s' in file '%s'" %
               (nodename, attributes, self.filename))
        return nodes[0] if nodes else None

    def get_nodes(self, nodename, attributes=None, root=None):
        if root is None:
            root = self.root
        nodes = []
        xpath = ".//"+nodename
        if attributes is not None:
            # xml.etree has limited support for xpath and does not allow more than
            # one attribute in an xpath query so we query seperately for each attribute
            # and create a result with the intersection of those lists
            for key, value in attributes.iteritems():
                xpath = ".//%s[@%s=\'%s\']" % (nodename, key, value)
                logger.debug("xpath is %s"%xpath)
                newnodes = root.findall(xpath)
                if not nodes:
                    nodes = newnodes
                else:
                    for node in nodes[:]:
                        if node not in newnodes:
                            nodes.remove(node)
                if not nodes:
                    return []
        else:
            nodes = root.findall(xpath)

        return nodes

    def add_child(self, node, root=None):
        """
        Add element node to self at root
        """
        if root is None:
            root = self.root
        self.root.append(node)

    def get_value(self, item, resolved=True):
        """
        get_value is expected to be defined by the derived classes, if you get here it is an error.
        """
        logger.debug("Get Value for "+item)
        result = None
        if item in self.lookups:
            result = self.lookups[item]

        if result is None:
            logger.debug("No value available for item '%s'" % item)
        elif resolved:
            result = self.get_resolved_value(result)

        return result

    def set_value(self, vid, value, ignore_type=True):
        """
        ignore_type is not used in this flavor
        """
        valnodes = self.get_nodes(vid)
        if valnodes:
            for node in valnodes:
                node.text = value
        else:
            print "DEBUG: adding %s %s to lookup" %(vid,value)
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
        logger.debug("raw_value %s" % raw_value)
        reference_re = re.compile(r'\$(\w+)')
        env_ref_re   = re.compile(r'\$ENV\{(\w+)\}')
        item_data = raw_value

        if item_data is None:
            return None

        if type(item_data) is not str:
            return item_data

        for m in env_ref_re.finditer(item_data):
            logger.debug("look for %s in env" % item_data)
            env_var = m.groups()[0]
            expect(env_var in os.environ, "Undefined env var '%s'" % env_var)
            item_data = item_data.replace(m.group(), os.environ[env_var])

        for m in reference_re.finditer(item_data):
            var = m.groups()[0]
            logger.debug("find: %s" % var)
            ref = self.get_value(var)
            if ref is not None:
                logger.debug("resolve: %s" % ref)
                item_data = item_data.replace(m.group(), str(self.get_resolved_value(ref)))
            elif var in os.environ:
                logging.warn("resolve from env: %s" % var)
                item_data = item_data.replace(m.group(), os.environ[var])

        return item_data
