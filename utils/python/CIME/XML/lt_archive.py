"""
Interface to the config_lt_archive.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.utils import expect
from CIME.XML.files import Files

logger = logging.getLogger(__name__)

class LTArchive(GenericXML):

    def __init__(self, machine, files=None, infile=None):
        """
        initialize an object
        """
        if (infile is None):
            if files is None:
                files = Files()
            infile = files.get_value("LTARCHIVE_SPEC_FILE", resolved=False)
            infile = files.get_resolved_value(infile)

        GenericXML.__init__(self, infile)

        self.machine = machine

    def get_lt_archive_args(self):
        """
        Return the lt archive args
        """
        nodes = self.get_nodes("machine")
        default = None
        for node in nodes:
            if node.get("name") == self.machine:
                return self.get_node("lt_archive_args", root=node).text
            elif node.get("name") == "default":
                default = self.get_node("lt_archive_args", root=node).text

        return default
