"""
Common interface to XML files, this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""
from CIME.XML.standard_module_setup import *
from CIME.utils import safe_copy

import xml.etree.ElementTree as ET
#pylint: disable=import-error
from distutils.spawn import find_executable
import getpass
import six
from copy import deepcopy
from collections import namedtuple

logger = logging.getLogger(__name__)

class _Element(object): # private class, don't want users constructing directly or calling methods on it

    def __init__(self, xml_element):
        self.xml_element = xml_element

    def __eq__(self, rhs):
        expect(isinstance(rhs, _Element), "Wrong type")
        return self.xml_element == rhs.xml_element # pylint: disable=protected-access

    def __ne__(self, rhs):
        expect(isinstance(rhs, _Element), "Wrong type")
        return self.xml_element != rhs.xml_element # pylint: disable=protected-access

    def __hash__(self):
        return hash(self.xml_element)

    def __deepcopy__(self, _):
        return _Element(deepcopy(self.xml_element))

class GenericXML(object):

    _FILEMAP = {}
    DISABLE_CACHING = False
    CacheEntry = namedtuple("CacheEntry", ["tree", "root", "modtime"])

    @classmethod
    def invalidate(cls, filename):
        if filename in cls._FILEMAP:
            del cls._FILEMAP[filename]

    def __init__(self, infile=None, schema=None, root_name_override=None, root_attrib_override=None, read_only=True):
        """
        Initialize an object
        """
        logger.debug("Initializing {}".format(infile))
        self.tree = None
        self.root = None
        self.locked = False
        self.read_only = read_only
        self.filename = infile
        self.needsrewrite = False
        if infile is None:
            return

        if os.path.isfile(infile) and os.access(infile, os.R_OK):
            # If file is defined and exists, read it
            self.read(infile, schema)
        else:
            # if file does not exist create a root xml element
            # and set it's id to file
            expect(not self.read_only, "Makes no sense to have empty read-only file")
            logger.debug("File {} does not exists.".format(infile))
            expect("$" not in infile,"File path not fully resolved {}".format(infile))

            root = _Element(ET.Element("xml"))

            if root_name_override:
                self.root = self.make_child(root_name_override, root=root, attributes=root_attrib_override)
            else:
                self.root = self.make_child("file", root=root, attributes={"id":os.path.basename(infile), "version":"2.0"})

            self.tree = ET.ElementTree(root)

            self._FILEMAP[infile] = self.CacheEntry(self.tree, self.root, 0.0)

    def read(self, infile, schema=None):
        """
        Read and parse an xml file into the object
        """
        cached_read = False
        if not self.DISABLE_CACHING and infile in self._FILEMAP:
            timestamp_cache = self._FILEMAP[infile].modtime
            timestamp_file  = os.path.getmtime(infile)
            if timestamp_file == timestamp_cache:
                logger.debug("read (cached): {}".format(infile))
                expect(self.read_only or not self.filename or not self.needsrewrite, "Reading into object marked for rewrite, file {}"
                       .format(self.filename))
                self.tree, self.root, _ = self._FILEMAP[infile]
                cached_read = True

        if not cached_read:
            logger.debug("read: {}".format(infile))
            file_open = (lambda x: open(x, 'r', encoding='utf-8')) if six.PY3 else (lambda x: open(x, 'r'))
            with file_open(infile) as fd:
                self.read_fd(fd)

            if schema is not None and self.get_version() > 1.0:
                self.validate_xml_file(infile, schema)

            logger.debug("File version is {}".format(str(self.get_version())))

            self._FILEMAP[infile] = self.CacheEntry(self.tree, self.root, os.path.getmtime(infile))

    def read_fd(self, fd):
        expect(self.read_only or not self.filename or not self.needsrewrite, "Reading into object marked for rewrite, file {}"               .format(self.filename))
        if self.tree:
            addroot = _Element(ET.parse(fd).getroot())
            read_only = self.read_only
            # we need to override the read_only mechanism here to append the xml object
            if read_only:
                self.read_only = False
            if addroot.xml_element.tag == self.name(self.root):
                for child in self.get_children(root=addroot):
                    self.add_child(child)
            else:
                self.add_child(addroot)
            if read_only:
                self.read_only = True
        else:
            self.tree = ET.parse(fd)
            self.root = _Element(self.tree.getroot())

    def lock(self):
        """
        A subclass is doing caching, we need to lock the tree structure
        in order to avoid invalidating cache.
        """
        self.locked = True

    def unlock(self):
        self.locked = False

    def change_file(self, newfile, copy=False):
        if copy:
            new_case = os.path.dirname(newfile)
            if not os.path.exists(new_case):
                os.makedirs(new_case)
            safe_copy(self.filename, newfile)

        self.tree = None
        self.filename = newfile
        self.read(newfile)

    #
    # API for individual node operations
    #

    def get(self, node, attrib_name, default=None):
        return node.xml_element.get(attrib_name, default=default)

    def has(self, node, attrib_name):
        return attrib_name in node.xml_element.attrib

    def set(self, node, attrib_name, value):
        expect(not self.read_only, "locked")
        if attrib_name == "id":
            expect(not self.locked, "locked")
        if self.get(node, attrib_name) != value:
            self.needsrewrite = True
            return node.xml_element.set(attrib_name, value)

    def pop(self, node, attrib_name):
        expect(not self.read_only, "locked")
        if attrib_name == "id":
            expect(not self.locked, "locked")
        self.needsrewrite = True
        return node.xml_element.attrib.pop(attrib_name)

    def attrib(self, node):
        # Return a COPY. We do not want clients making changes directly
        return None if node.xml_element.attrib is None else dict(node.xml_element.attrib)

    def set_name(self, node, name):
        expect(not self.read_only, "locked")
        if node.xml_element.tag != name:
            self.needsrewrite = True
            node.xml_element.tag = name

    def set_text(self, node, text):
        expect(not self.read_only, "locked")
        if node.xml_element.text != text:
            node.xml_element.text = text
            self.needsrewrite = True

    def name(self, node):
        return node.xml_element.tag

    def text(self, node):
        return node.xml_element.text

    def add_child(self, node, root=None, position=None):
        """
        Add element node to self at root
        """
        expect(not self.locked and not self.read_only, "locked")
        self.needsrewrite = True
        root = root if root is not None else self.root
        if position is not None:
            root.xml_element.insert(position, node.xml_element)
        else:
            root.xml_element.append(node.xml_element)

    def copy(self, node):
        return deepcopy(node)

    def remove_child(self, node, root=None):
        expect(not self.locked and not self.read_only, "locked")
        self.needsrewrite = True
        root = root if root is not None else self.root
        root.xml_element.remove(node.xml_element)

    def make_child(self, name, attributes=None, root=None, text=None):
        expect(not self.locked and not self.read_only, "locked")
        root = root if root is not None else self.root
        self.needsrewrite = True
        if attributes is None:
            node = _Element(ET.SubElement(root.xml_element, name))
        else:
            node = _Element(ET.SubElement(root.xml_element, name, attrib=attributes))

        if text:
            self.set_text(node, text)

        return node

    def get_children(self, name=None, attributes=None, root=None):
        """
        This is the critical function, its interface and performance are crucial.

        You can specify attributes={key:None} if you want to select chilren
        with the key attribute but you don't care what its value is.
        """
        root = root if root is not None else self.root
        children = []
        for child in root.xml_element:
            if name is not None:
                if child.tag != name:
                    continue

            if attributes is not None:
                if child.attrib is None:
                    continue
                else:
                    match = True
                    for key, value in attributes.items():
                        if key not in child.attrib:
                            match = False
                            break
                        elif value is not None:
                            if child.attrib[key] != value:
                                match = False
                                break

                    if not match:
                        continue

            children.append(_Element(child))

        return children

    def get_child(self, name=None, attributes=None, root=None, err_msg=None):
        children = self.get_children(root=root, name=name, attributes=attributes)
        expect(len(children) == 1, err_msg if err_msg else "Expected one child")
        return children[0]

    def get_optional_child(self, name=None, attributes=None, root=None, err_msg=None):
        children = self.get_children(root=root, name=name, attributes=attributes)
        expect(len(children) <= 1, err_msg if err_msg else "Multiple matches")
        return children[0] if children else None

    def get_element_text(self, element_name, attributes=None, root=None):
        element_node = self.get_optional_child(name=element_name, attributes=attributes, root=root)
        if element_node is not None:
            return self.text(element_node)
        return None

    def set_element_text(self, element_name, new_text, attributes=None, root=None):
        element_node = self.get_optional_child(name=element_name, attributes=attributes, root=root)
        if element_node is not None:
            self.set_text(element_node, new_text)
            return new_text
        return None

    def to_string(self, node, method="xml", encoding="us-ascii"):
        return ET.tostring(node, method=method, encoding=encoding)

    #
    # API for operations over the entire file
    #

    def get_version(self):
        version = self.get(self.root, "version")
        version = 1.0 if version is None else float(version)
        return version

    def write(self, outfile=None, force_write=False):
        """
        Write an xml file from data in self
        """
        timestamp_cache = self._FILEMAP[self.filename].modtime
        if timestamp_cache != 0.0:
            timestamp_file  = os.path.getmtime(self.filename)
            expect(timestamp_file == timestamp_cache,
                   "File {} appears to have changed without a corresponding invalidation, modtimes {:0.2f} != {:0.2f}".format(self.filename, timestamp_cache, timestamp_file))

        if not (self.needsrewrite or force_write):
            return

        if outfile is None:
            outfile = self.filename

        logger.debug("write: " + (outfile if isinstance(outfile, six.string_types) else str(outfile)))

        xmlstr = self.get_raw_record()

        # xmllint provides a better format option for the output file
        xmllint = find_executable("xmllint")
        if xmllint is not None:
            if isinstance(outfile, six.string_types):
                run_cmd_no_fail("{} --format --output {} -".format(xmllint, outfile), input_str=xmlstr)
            else:
                outfile.write(run_cmd_no_fail("{} --format -".format(xmllint), input_str=xmlstr))

        else:
            with open(outfile,'w') as xmlout:
                xmlout.write(xmlstr)

        self._FILEMAP[self.filename] = self.CacheEntry(self.tree, self.root, os.path.getmtime(self.filename))

        self.needsrewrite = False

    def scan_child(self, nodename, attributes=None, root=None):
        """
        Get an xml element matching nodename with optional attributes.

        Error unless exactly one match.
        """

        nodes = self.scan_children(nodename, attributes=attributes, root=root)

        expect(len(nodes) == 1, "Incorrect number of matches, {:d}, for nodename '{}' and attrs '{}' in file '{}'".format(len(nodes), nodename, attributes, self.filename))
        return nodes[0]

    def scan_optional_child(self, nodename, attributes=None, root=None):
        """
        Get an xml element matching nodename with optional attributes.

        Return None if no match.
        """
        nodes = self.scan_children(nodename, attributes=attributes, root=root)

        expect(len(nodes) <= 1, "Multiple matches for nodename '{}' and attrs '{}' in file '{}'".format(nodename, attributes, self.filename))
        return nodes[0] if nodes else None

    def scan_children(self, nodename, attributes=None, root=None):

        logger.debug("(get_nodes) Input values: {}, {}, {}, {}".format(self.__class__.__name__, nodename, attributes, root))

        if root is None:
            root = self.root
        nodes = []

        xpath = ".//" + (nodename if nodename else "")

        if attributes:
            # xml.etree has limited support for xpath and does not allow more than
            # one attribute in an xpath query so we query seperately for each attribute
            # and create a result with the intersection of those lists

            for key, value in attributes.items():
                if value is None:
                    xpath = ".//{}[@{}]".format(nodename, key)
                else:
                    xpath = ".//{}[@{}=\'{}\']".format(nodename, key, value)

                logger.debug("xpath is {}".format(xpath))

                try:
                    newnodes = root.xml_element.findall(xpath)
                except Exception as e:
                    expect(False, "Bad xpath search term '{}', error: {}".format(xpath, e))

                if not nodes:
                    nodes = newnodes
                else:
                    for node in nodes[:]:
                        if node not in newnodes:
                            nodes.remove(node)
                if not nodes:
                    return []

        else:
            logger.debug("xpath: {}".format(xpath))
            nodes = root.xml_element.findall(xpath)

        logger.debug("Returning {} nodes ({})".format(len(nodes), nodes))

        return [_Element(node) for node in nodes]

    def get_value(self, item, attribute=None, resolved=True, subgroup=None): # pylint: disable=unused-argument
        """
        get_value is expected to be defined by the derived classes, if you get here
        the value was not found in the class.
        """
        logger.debug("Get Value for " + item)
        return None

    def get_values(self, vid, attribute=None, resolved=True, subgroup=None):# pylint: disable=unused-argument
        logger.debug("Get Values for " + vid)
        return []

    def set_value(self, vid, value, subgroup=None, ignore_type=True): # pylint: disable=unused-argument
        """
        ignore_type is not used in this flavor
        """
        valnodes = self.get_children(vid)
        for node in valnodes:
            self.set_text(node, value)

        return value if valnodes else None

    def get_resolved_value(self, raw_value, allow_unresolved_envvars=False):
        """
        A value in the xml file may contain references to other xml
        variables or to environment variables. These are refered to in
        the perl style with $name and $ENV{name}.

        >>> obj = GenericXML()
        >>> os.environ["FOO"] = "BAR"
        >>> os.environ["BAZ"] = "BARF"
        >>> obj.get_resolved_value("one $ENV{FOO} two $ENV{BAZ} three")
        'one BAR two BARF three'
        >>> obj.get_resolved_value("2 + 3 - 1")
        '4'
        >>> obj.get_resolved_value("0001-01-01")
        '0001-01-01'
        >>> obj.get_resolved_value("$SHELL{echo hi}") == 'hi'
        True
        """
        logger.debug("raw_value {}".format(raw_value))
        reference_re = re.compile(r'\${?(\w+)}?')
        env_ref_re   = re.compile(r'\$ENV\{(\w+)\}')
        shell_ref_re = re.compile(r'\$SHELL\{([^}]+)\}')
        math_re = re.compile(r'\s[+-/*]\s')
        item_data = raw_value

        if item_data is None:
            return None

        if not isinstance(item_data, six.string_types):
            return item_data

        for m in env_ref_re.finditer(item_data):
            logger.debug("look for {} in env".format(item_data))
            env_var = m.groups()[0]
            env_var_exists = env_var in os.environ
            if not allow_unresolved_envvars:
                expect(env_var_exists, "Undefined env var '{}'".format(env_var))
            if env_var_exists:
                item_data = item_data.replace(m.group(), os.environ[env_var])

        for s in shell_ref_re.finditer(item_data):
            logger.debug("execute {} in shell".format(item_data))
            shell_cmd = s.groups()[0]
            item_data = item_data.replace(s.group(), run_cmd_no_fail(shell_cmd))

        for m in reference_re.finditer(item_data):
            var = m.groups()[0]
            logger.debug("find: {}".format(var))
            # The overridden versions of this method do not simply return None
            # so the pylint should not be flagging this
            ref = self.get_value(var) # pylint: disable=assignment-from-none

            if ref is not None:
                logger.debug("resolve: " + str(ref))
                item_data = item_data.replace(m.group(), self.get_resolved_value(str(ref)))
            elif var == "CIMEROOT":
                cimeroot = get_cime_root()
                item_data = item_data.replace(m.group(), cimeroot)
            elif var == "SRCROOT":
                srcroot = os.path.join(get_cime_root(),"..")
                item_data = item_data.replace(m.group(), srcroot)
            elif var == "USER":
                item_data = item_data.replace(m.group(), getpass.getuser())

        if math_re.search(item_data):
            try:
                tmp = eval(item_data)
            except:
                tmp = item_data
            item_data = str(tmp)

        return item_data

    def validate_xml_file(self, filename, schema):
        """
        validate an XML file against a provided schema file using pylint
        """
        expect(os.path.isfile(filename),"xml file not found {}".format(filename))
        expect(os.path.isfile(schema),"schema file not found {}".format(schema))
        xmllint = find_executable("xmllint")
        if xmllint is not None:
            logger.debug("Checking file {} against schema {}".format(filename, schema))
            run_cmd_no_fail("{} --noout --schema {} {}".format(xmllint, schema, filename))
        else:
            logger.warning("xmllint not found, could not validate file {}".format(filename))

    def get_raw_record(self, root=None):
        logger.debug("writing file {}".format(self.filename))
        if root is None:
            root = self.root
        try:
            xmlstr = ET.tostring(root.xml_element)
        except ET.ParseError as e:
            ET.dump(root.xml_element)
            expect(False, "Could not write file {}, xml formatting error '{}'".format(self.filename, e))
        return xmlstr

    def get_id(self):
        xmlid = self.get(self.root, "id")
        if xmlid is not None:
            return xmlid
        return self.name(self.root)
