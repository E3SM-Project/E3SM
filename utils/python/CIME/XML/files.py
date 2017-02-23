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
        '$CIMEROOT/config/config_headers.xml'
        """
        cimeroot = get_cime_root()
        infile = os.path.join(cimeroot, "config", get_model(), "config_files.xml")
        schema = os.path.join(cimeroot, "config", "xml_schemas", "entry_id.xsd")
        EntryID.__init__(self, infile, schema=schema)

    def get_schema(self, nodename):
        node = self.get_optional_node("entry", {"id":nodename})
        schemanode = self.get_optional_node("schema", root=node)
        if schemanode is not None:
            logger.debug("Found schema for %s"%(nodename))
            return self.get_resolved_value(schemanode.text)
        return None

    def get_components(self, nodename):
        node = self.get_optional_node("entry", {"id":nodename})
        valnodes = self.get_nodes("value", root=node)
        values = []
        for valnode in valnodes:
            value = valnode.get("component")
            values.append(value)
        return values
