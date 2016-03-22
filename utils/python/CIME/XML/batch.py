"""
Interface to the config_batch.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Batch(GenericXML):

    def __init__(self, infile=None, batch_system=None):
        """
        initialize an object
        """
        if infile is None:
            infile = os.path.join(get_cime_root(), "cime_config", get_model(), "machines", "config_batch.xml")

        GenericXML.__init__(self, infile)

        self.batch_system = None
        self.name         = batch_system

    def get_batch_system(self):
        """
        Return the name of the batch system
        """
        return self.name

    def get_node(self, nodename, attributes=None):
        """
        Return data on a node for a batch system
        """
        expect(self.batch_system is not None, "Batch system not set, use parent get_node?")
        return GenericXML.get_optional_node(self, nodename, attributes, root=self.batch_system)

    def get_optional_node(self, nodename, attributes=None):
        """
        Return data on a node for a batch system
        """
        expect(self.batch_system is not None, "Batch system not set, use parent get_node?")
        return GenericXML.get_optional_node(self, nodename, attributes, root=self.batch_system)

    def set_batch_system(self, batch_system):
        """
        Sets the batch system block in the Batch object
        """
        if self.name != batch_system:
            self.batch_system = GenericXML.get_optional_node(self, "batch_system", {"type" : batch_system})
            expect(self.batch_system is not None, "No batch system '%s' found" % batch_system)
            self.name = batch_system

        return batch_system

    def get_value(self, name, resolved=True):
        """
        Get Value of fields in the config_batch.xml file
        """
        expect(self.batch_system is not None, "Batch object has no batch system defined")
        value = None

        node = self.get_optional_node(name)
        if node is not None:
            value = node.text

        if value is None:
            # if all else fails
            value = GenericXML.get_value(self, name)

        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

        return value
