"""
Interface to the config_component.xml files.  This class inherits from EntryID.py
"""
from standard_module_setup import *

from entry_id import EntryID
from CIME.utils import expect, get_cime_root, get_model
from CIME.XML.files import Files

logger = logging.getLogger(__name__)

class Component(EntryID):

    def __init__(self, infile=None):
        """
        initialize an object
        """
        if infile is None:
            files = Files()
            infile = files.get_value("CONFIG_DRV_FILE")

        EntryID.__init__(self,infile)

    def get_value(self, name, attribute={}, resolved=False):
        return EntryID.get_value(self, name, attribute, resolved)

    def get_valid_model_components(self):
        components = []
        compnode = self.get_node("components")
        comps = self.get_nodes("comp", root=compnode)
        for comp in comps:
            components.append(comp.text)
        return components
