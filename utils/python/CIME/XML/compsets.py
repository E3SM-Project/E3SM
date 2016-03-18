"""
Common interface to XML files which follow the compsets format,
"""

from standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type
from generic_xml import GenericXML

logger = logging.getLogger(__name__)

class Compsets(GenericXML):

    def __init__(self, infile=None):
        GenericXML.__init__(self, infile)
        self.groups={}

    def get_value(self, name, attribute={}, resolved=False):
        if name == "help":
            rootnode = self.get_node("help")
            helptext = rootnode.text
            return helptext
        else:
            compsets = {}
            nodes = self.get_nodes(nodename="compset")
            for node in nodes:
                for child in node:
                    logger.debug ("Here child is %s with value %s"%(child.tag,child.text))
                    if child.tag == "alias": 
                        alias = child.text
                    if child.tag == "lname":
                        lname = child.text
                compsets[alias] = lname
            return compsets

    def print_values(self):
        help_text = self.get_value(name="help")
        compsets_text = self.get_value("names")
        print help_text

        print "{:<30} {:<35} ".format('Compset Alias','Compset Longname')
        print "{:<30} {:<35} ".format('=============','================')
        for v in compsets_text.iteritems():
            label, definition = v
            print "{:<30} {:<35} ".format(label, definition)

