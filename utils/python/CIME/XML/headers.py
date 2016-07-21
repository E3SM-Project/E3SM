"""
Interface to the config_headers.xml file.  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.entry_id import EntryID
from CIME.XML.files import Files

logger = logging.getLogger(__name__)

class Headers(EntryID):
    def __init__(self,infile=None):
        """
        initialize an object

        >>> files = Files()
        >>> files.get_value('CASEFILE_HEADERS',resolved=False)
        '$CIMEROOT/cime_config/config_headers.xml'
        """
        if infile is None:
            files = Files()
            infile = files.get_value('CASEFILE_HEADERS', resolved=True)
        EntryID.__init__(self, infile)

    def get_header_node(self, fname):
        fnode = self.get_node("file", attributes={"name" : fname})
        headernode = self.get_node("header", root=fnode)
        return headernode
