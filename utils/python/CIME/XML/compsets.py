"""
Common interface to XML files which follow the compsets format,
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type
from CIME.XML.generic_xml import GenericXML

logger = logging.getLogger(__name__)

class Compsets(GenericXML):

    def __init__(self, infile=None):
        GenericXML.__init__(self, infile)
        self.groups={}

    def get_compset_match(self, name):
        nodes = self.get_nodes("compset")
        for node in nodes:
            alias = self.get_node("alias",root=node)
            lname = self.get_node("lname",root=node)
            if alias.text == name or lname.text == name:
                logger.debug("Found node match with alias: %s and lname: %s" % (alias.text, lname.text))
                return lname.text

    def get_value(self, name, attribute={}, resolved=False, subgroup=None):
        expect(subgroup is None, "This class does not support subgroups")
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
        logger.info(" %s " %help_text)

        logger.info("       --------------------------------------")
        logger.info("       Compset Short Name: Compset Long Name ")
        logger.info("       --------------------------------------")
        for v in compsets_text.iteritems():
            label, definition = v
            logger.info("   %20s : %s" %(label, definition))

