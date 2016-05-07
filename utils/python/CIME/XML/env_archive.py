"""
Interface to the env_archive.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.generic_xml import GenericXML
from CIME.XML.archive import Archive
from CIME.XML.headers import Headers

logger = logging.getLogger(__name__)

class EnvArchive(GenericXML):

    def __init__(self, case_root=None, infile="env_archive.xml"):
        """
        initialize an object interface to file env_archive.xml in the case directory
        """
        logger.debug("Case_root = %s" , case_root)

        # Check/Build path to env_archive.xml
        if case_root is None:
            case_root = os.getcwd()

        if os.path.isabs(infile):
            fullpath = infile
        else:
            fullpath = os.path.join(case_root, infile)

        # Initialize self
        # If env_archive.xml file does not exists in case directory read default from config
        GenericXML.__init__(self, fullpath)

        # The following creates the CASEROOT/env_archive.xml contents in self.root
        if not os.path.isfile(fullpath):
            headerobj = Headers()
            headernode = headerobj.get_header_node(os.path.basename(fullpath))
            self.root.append(headernode)
            archive = Archive()
            self.root.append(archive.root)

    def get_values(self, item, attribute={}, resolved=True, subgroup=None):
        """Returns the value as a string of the first xml element with item as attribute value.
        <elememt_name attribute='attribute_value>value</element_name>"""

        logger.debug("(get_values) Input values: %s , %s , %s , %s , %s" , self.__class__.__name__, item, attribute, resolved, subgroup)

        nodes   = [] # List of identified xml elements
        results = [] # List of identified parameters

        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        if item :
            logger.debug("Searching for %s" , item )
            nodes = self.get_nodes("*",{"name" : item})
        else :
            # Return all nodes
            logger.debug("Retrieving all parameter")
            nodes = self.get_nodes("comp_archive_spec")

        logger.debug("Number found %s" , nodes.__len__())
        # Return value for first occurence of node with attribute value = item
        for node in nodes:
            logger.debug("Checking node %s" , node.tag )
            group   = node.attrib['name']
            rootdir = node.find('./rootdir').text

            for file_ in node.findall('./files/file_extension') :

                keep    = file_.find('./keep_last_in_rundir')

                val     = keep.text
                attr    = rootdir + "/" + file_.find('./subdir').text + "/" + file_.get('regex_suffix')
                t       = self._get_type(node)
                desc    = self._get_description(keep)
                #default = super(EnvBase , self).get_default(node)
                default = self._get_default(keep)
                filename = self.filename

                v = { 'group' : group , 'attribute' : attr , 'value' : val , 'type' : t , 'description' : desc , 'default' : default , 'file' : filename}
                results.append(v)

        return results

    def _get_type(self, node=None) :
        return None

    def _get_description(self, node=None):

        description = None

        if node.tag == "keep_last_in_rundir" :
            if node.get('description') :
                description = node.get('description')
            elif node.find('./description') :
                description = node.find('./description').text
            else:
                description = "Keep last in rundir"

        return description

    def _get_default(self, node=None):
        return None

    def _get_archive():
        pass
