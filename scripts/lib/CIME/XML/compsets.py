"""
Common interface to XML files which follow the compsets format,
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.entry_id import EntryID
from CIME.XML.files import Files


logger = logging.getLogger(__name__)

class Compsets(GenericXML):

    def __init__(self, infile=None, files=None):
        if files is None:
            files = Files()
        schema = files.get_schema("COMPSETS_SPEC_FILE")
        GenericXML.__init__(self, infile, schema=schema)
        self.groups={}

    def get_compset_match(self, name):
        """
        science support is used in cesm to determine if this compset and grid
        is scientifically supported.   science_support is returned as an array of grids for this compset
        """
        nodes = self.get_nodes("compset")
        alias = None
        lname = None
        science_support = []
        for node in nodes:
            alias = self.get_element_text("alias",root=node)
            lname = self.get_element_text("lname",root=node)
            if alias == name or lname == name:
                science_support_nodes = self.get_nodes("science_support", root=node)
                for node in science_support_nodes:
                    science_support.append(node.get("grid"))

                logger.debug("Found node match with alias: %s and lname: %s" % (alias, lname))
                return (lname, alias, science_support)
        return (None, None, False)

    def get_compset_var_settings(self, compset, grid):
        '''
        Variables can be set in config_compsets.xml in entry id settings with compset and grid attributes
        find and return id value pairs here
        '''
        nodes = self.get_nodes("entry")
        # Get an empty entryid obj to use
        entryidobj = EntryID()
        result = []
        for node in nodes:
            value = entryidobj.get_default_value(node, {"grid":grid, "compset":compset})
            if value is not None:
                result.append((node.get("id"), value))
        return result

    def get_value(self, name, attribute=None, resolved=False, subgroup=None):
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

    def return_all_values(self):
        all_compsets = dict()
        help_text = self.get_value(name="help")
        compsets_text = self.get_value("names")
        for v in compsets_text.iteritems():
            label, definition = v
            all_compsets[label] = (help_text, definition)

        return all_compsets

