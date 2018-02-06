"""
Interface to the env_archive.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvArchive(EnvBase):

    def __init__(self, case_root=None, infile="env_archive.xml"):
        """
        initialize an object interface to file env_archive.xml in the case directory
        """
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_archive.xsd")
        EnvBase.__init__(self, case_root, infile, schema=schema)

    def get_entries(self):
        return self.get_children('comp_archive_spec')

    def get_entry(self, compname):
        components = self.get_optional_child('components')
        return None if components is None else self.get_optional_child('comp_archive_spec', attributes={"compname":compname}, root=components)

    def get_entry_info(self, archive_entry):
        compname = self.get(archive_entry, 'compname')
        compclass = self.get(archive_entry, 'compclass')
        return compname,compclass

    def get_entry_value(self, name, archive_entry):
        node = self.get_optional_child(name, root=archive_entry)
        if node is not None:
            return self.text(node)
        return None

    def get_rest_file_extensions(self, archive_entry):
        file_extensions = []
        nodes = self.get_children('rest_file_extension', root=archive_entry)
        for node in nodes:
            file_extensions.append(self.text(node))
        return file_extensions

    def get_hist_file_extensions(self, archive_entry):
        file_extensions = []
        nodes = self.get_children('hist_file_extension', root=archive_entry)
        for node in nodes:
            file_extensions.append(self.text(node))
        return file_extensions

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
