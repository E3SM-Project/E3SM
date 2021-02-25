"""
Interface to the env_archive.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.archive_base import ArchiveBase
from CIME.XML.env_base import EnvBase

logger = logging.getLogger(__name__)
# pylint: disable=super-init-not-called
class EnvArchive(ArchiveBase,EnvBase):

    def __init__(self, case_root=None, infile="env_archive.xml", read_only=False):
        """
        initialize an object interface to file env_archive.xml in the case directory
        """
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_archive.xsd")
        EnvBase.__init__(self, case_root, infile, schema=schema, read_only=read_only)


    def get_entries(self):
        return self.get_children('comp_archive_spec')

    def get_entry_info(self, archive_entry):
        compname = self.get(archive_entry, 'compname')
        compclass = self.get(archive_entry, 'compclass')
        return compname,compclass

    def get_rpointer_contents(self, archive_entry):
        rpointer_items = []
        rpointer_nodes = self.get_children('rpointer', root=archive_entry)
        for rpointer_node in rpointer_nodes:
            file_node = self.get_child('rpointer_file', root=rpointer_node)
            content_node = self.get_child('rpointer_content', root=rpointer_node)
            rpointer_items.append([self.text(file_node),self.text(content_node)])
        return rpointer_items

    def get_type_info(self, vid):
        return "char"
