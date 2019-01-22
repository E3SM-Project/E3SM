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

    def get_compset_match(self, name):
        """
        science support is used in cesm to determine if this compset and grid
        is scientifically supported.   science_support is returned as an array of grids for this compset
        """
        nodes = self.get_children("compset")
        alias = None
        lname = None

        science_support = []

        for node in nodes:
            alias = self.get_element_text("alias",root=node)
            lname = self.get_element_text("lname",root=node)
            if alias == name or lname == name:
                science_support_nodes = self.get_children("science_support", root=node)
                for snode in science_support_nodes:
                    science_support.append(self.get(snode, "grid"))
                logger.debug("Found node match with alias: {} and lname: {}".format(alias, lname))
                return (lname, alias, science_support)
        return (None, None, [False])

    def get_compset_var_settings(self, compset, grid):
        '''
        Variables can be set in config_compsets.xml in entry id settings with compset and grid attributes
        find and return id value pairs here
        '''
        entries = self.get_optional_child("entries")
        result = []
        if entries is not None:
            nodes = self.get_children("entry", root=entries)
            # Get an empty entryid obj to use
            entryidobj = EntryID()
            for node in nodes:
                value = entryidobj.get_default_value(node, {"grid":grid, "compset":compset})
                if value is not None:
                    result.append((self.get(node, "id"), value))

        return result

    #pylint: disable=arguments-differ
    def get_value(self, name, attribute=None, resolved=False, subgroup=None):
        expect(subgroup is None, "This class does not support subgroups")
        if name == "help":
            rootnode = self.get_child("help")
            helptext = self.text(rootnode)
            return helptext
        else:
            compsets = {}
            nodes = self.get_children("compset")
            for node in nodes:
                for child in node:
                    logger.debug ("Here child is {} with value {}".format(self.name(child),self.text(child)))
                    if self.name(child) == "alias":
                        alias = self.text(child)
                    if self.name(child) == "lname":
                        lname = self.text(child)
                compsets[alias] = lname
            return compsets

    def print_values(self, arg_help=True):
        help_text = self.get_value(name="help")
        compsets = self.get_children("compset")
        if arg_help:
            logger.info(" {} ".format(help_text))

        logger.info("       --------------------------------------")
        logger.info("       Compset Alias: Compset Long Name ")
        logger.info("       --------------------------------------")
        for compset in compsets:
            logger.info("   {:20} : {}".format(self.text(self.get_child("alias",root=compset)),
                                               self.text(self.get_child("lname", root=compset))))

    def get_compset_longnames(self):
        compset_nodes = self.get_children("compset")
        longnames = []
        for comp in compset_nodes:
            longnames.append(self.text(self.get_child("lname", root=comp)))
        return(longnames)
