"""
Interface to the config_inputdata.xml file.  This class inherits from GenericXML.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Inputdata(GenericXML):

    def __init__(self, infile=None, files=None):
        """
        initialize a files object given input pes specification file
        """
        if files is None:
            files = Files()
        if infile is None:
            infile = files.get_value("INPUTDATA_SPEC_FILE")
        schema = files.get_schema("INPUTDATA_SPEC_FILE")
        logger.debug("DEBUG: infile is {}".format(infile))
        GenericXML.__init__(self, infile, schema=schema)

        self._servernode = None

    def get_next_server(self):
        protocol = None
        address = None
        servernodes = self.get_children("server")
        if self._servernode is None:
            self._servernode = servernodes[0]
        else:
            prevserver = self._servernode
            for i, node in enumerate(servernodes):
                if self._servernode == node and len(servernodes)>i+1:
                    self._servernode = servernodes[i+1]
                    break
            if prevserver is not None and self._servernode == prevserver:
                self._servernode = None

        if self._servernode is not None:
            protocol = self.text(self.get_child("protocol", root = self._servernode))
            address =  self.text(self.get_child("address", root = self._servernode))

        return protocol, address
