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

    def get_file_node_text(self, attnames, archive_entry):
        nodes = []
        textvals = []
        for attname in attnames:
            nodes.extend(self.get_children(attname, root=archive_entry))
        for node in nodes:
            textvals.append(self.text(node))
        return textvals

    def get_rest_file_extensions(self, archive_entry):
        return self.get_file_node_text(['rest_file_extension'],archive_entry)

    def get_rest_file_regex(self, archive_entry):
        return self.get_file_node_text(['rest_file_regex'],archive_entry)

    def get_hist_file_extensions(self, archive_entry):
        return self.get_file_node_text(['hist_file_extension'],archive_entry)

    def get_hist_file_regex(self, archive_entry):
        return self.get_file_node_text(['hist_file_regex'],archive_entry)

    def get_entry_value(self, name, archive_entry):
        node = self.get_optional_child(name, root=archive_entry)
        if node is not None:
            return self.text(node)
        return None
