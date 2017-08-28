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
                for snode in science_support_nodes:
                    science_support.append(snode.get("grid"))
                logger.debug("Found node match with alias: {} and lname: {}".format(alias, lname))
                return (lname, alias, science_support)
        return (None, None, [False])

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

    #pylint: disable=arguments-differ
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
                    logger.debug ("Here child is {} with value {}".format(child.tag,child.text))
                    if child.tag == "alias":
                        alias = child.text
                    if child.tag == "lname":
                        lname = child.text
                compsets[alias] = lname
            return compsets

    def print_values(self, arg_help=True):
        help_text = self.get_value(name="help")
        compsets_text = self.get_value("names")
        if arg_help:
            logger.info(" {} ".format(help_text))

        logger.info("       --------------------------------------")
        logger.info("       Compset Alias: Compset Long Name ")
        logger.info("       --------------------------------------")
        for key in sorted(compsets_text.iterkeys()):
            logger.info("   {:20} : {}".format(key, compsets_text[key]))

    def return_all_values(self):
        all_compsets = dict()
        science_compsets = dict()
        help_text = self.get_value(name="help")
        compsets_text = self.get_value("names")
        for key in sorted(compsets_text.iterkeys()):
            all_compsets[key] = compsets_text[key]

        # get the matching science support grids
        for alias in all_compsets.iterkeys():
            science_compsets[alias] = self.get_compset_match(alias)
            
        return help_text, all_compsets

