"""
Interface to the streams.xml style files.  This class inherits from GenericXML.py

stream files predate cime and so do not conform to entry id format
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Stream(GenericXML):

    def __init__(self, infile=None, files=None):
        """
        initialize an object
        """
        if files is None:
            files = Files()
        schema = None
        GenericXML.__init__(self, infile, schema=schema)

    def get_value(self, name, attribute=None, resolved=True, subgroup=None):
        """
        Get Value of fields in a stream.xml file
        """
        expect(subgroup is None, "This class does not support subgroups")
        value = None
        node = None
        if name.startswith("domain"):
            node = self.scan_child("domainInfo")
        elif name.startswith("data"):
            node = self.scan_child("fieldInfo")
        if name.endswith("filepath"):
            node = self.scan_child("filePath", root=node)
        elif name.endswith("filenames"):
            node = self.scan_child("fileNames", root=node)
        if node is not None:
            value = self.text(node).strip()

        if value is None:
            # if all else fails
            #pylint: disable=assignment-from-none
            value = GenericXML.get_value(self, name, attribute, resolved, subgroup)

        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

        return value
