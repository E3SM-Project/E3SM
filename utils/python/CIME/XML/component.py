"""
Interface to the config_component.xml files.  This class inherits from EntryID.py
"""
from standard_module_setup import *

from entry_id import EntryID
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Component(EntryID):

    def __init__(self,infile):
        """
        initialize an object
        """
        EntryID.__init__(self,infile)

    def get_value(self, name, attribute={}, resolved=False):
        if name == "component":
            components = []
            compnode = self.get_node("components")
            comps = self.get_nodes("comp", root=compnode)
            for comp in comps:
                components.append(comp.text)
            return components
        else:
            return EntryID.get_value(self, name, attribute, resolved)
