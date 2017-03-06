"""
Interface to the env_archive.xml file.  This class inherits from GenericXML
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

    def get_entries(self):
        return self.get_nodes('comp_archive_spec')

    def get_entry_info(self, archive_entry):
        compname = archive_entry.attrib['compname']
        compclass = archive_entry.attrib['compclass']
        return compname,compclass

    def get_entry_value(self, name, archive_entry):
        node = self.get_optional_node(name, root=archive_entry)
        return node.text

    def get_rest_file_extensions(self, archive_entry):
        file_extensions = []
        nodes = self.get_nodes('rest_file_extension', root=archive_entry)
        for node in nodes:
            file_extensions.append(node.text)
        return file_extensions

    def get_hist_file_extensions(self, archive_entry):
        file_extensions = []
        nodes = self.get_nodes('hist_file_extension', root=archive_entry)
        for node in nodes:
            file_extensions.append(node.text)
        return file_extensions

    def get_rpointer_contents(self, archive_entry):
        rpointer_items = []
        rpointer_nodes = self.get_nodes('rpointer', root=archive_entry)
        for rpointer_node in rpointer_nodes:
            file_node = self.get_node('rpointer_file', root=rpointer_node)
            content_node = self.get_node('rpointer_content', root=rpointer_node)
            rpointer_items.append([file_node.text,content_node.text])
        return rpointer_items
