"""
Base class for archive files.  This class inherits from generic_xml.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML

logger = logging.getLogger(__name__)

class ArchiveBase(GenericXML):

    def get_entry(self, compname):
        return self.scan_optional_child('comp_archive_spec',
                                        attributes={"compname":compname})


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

    def get_entry_value(self, name, archive_entry):
        node = self.get_optional_child(name, root=archive_entry)
        if node is not None:
            return self.text(node)
        return None
