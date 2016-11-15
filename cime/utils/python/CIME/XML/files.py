"""
Interface to the config_files.xml file.  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.entry_id import EntryID
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Files(EntryID):

    def __init__(self):
        """
        initialize an object

        >>> files = Files()
        >>> files.get_value('CASEFILE_HEADERS',resolved=False)
        '$CIMEROOT/cime_config/config_headers.xml'
        """
        infile = os.path.join(get_cime_root(), "cime_config", get_model(), "config_files.xml")
        EntryID.__init__(self, infile)

    def get_components(self, nodename):
        node = self.get_optional_node("entry", {"id":nodename})
        valnodes = self.get_nodes("value", root=node)
        values = []
        for valnode in valnodes:
            value = valnode.get("component")
            values.append(value)
        return values
