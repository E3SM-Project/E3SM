"""
Common interface to XML files which follow the grids format,
This is not an abstract class - but inherits from the abstact class GenericXML
"""

from standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type
from generic_xml import GenericXML

logger = logging.getLogger(__name__)

class Grids(GenericXML):

    def __init__(self, infile=None):
        GenericXML.__init__(self, infile)
        self.groups={}

    def _get_component_grids(self, name):
        gridRE = re.compile(r"[_]{0,1}[a-z]{1,2}%")
        component_grids = gridRE.split(name.text)[1:]
        return component_grids

    def print_values(self, long_output=None):
        # write out help message
        help_node = self.get_node("help")
        helptext = help_node.text
        logger.info("%s " %helptext)

        # if long output mode is on, then also obtain nodes for domains and maps
        if long_output is not None:
            domain_nodes = self.get_nodes(nodename="domain")
            gridmap_nodes = self.get_nodes(nodename="gridmap")

        # write out grid elements
        grid_nodes = self.get_nodes(nodename="grid")
        for grid_node in grid_nodes:
            lname   = grid_node.find("lname")
            sname   = grid_node.find("sname")
            support = grid_node.find("support")
            alias   = grid_node.find("alias")
            logger.info("------------------------------------------------")
            logger.info("model grid: %s" %lname.text)
            logger.info("------------------------------------------------")
            if sname is not None:
                logger.info("   short name: %s" %sname.text)
            if alias is not None:
                logger.info("   alias: %s" %alias.text)
            for attr, value in grid_node.items():
                if  attr == 'compset':
                    logger.info("   compset match: %s" % value)

            # in long_output mode add domain description and mapping fiels
            if long_output is not None:

                # get component grids (will contain duplicates)
                component_grids = self._get_component_grids(lname)

                # write out out only unique non-null component grids
                for domain in list(set(component_grids)):
                    if domain != 'null':
                        logger.info("   domain is %s" %domain)
                        domain_node = self.get_node(nodename="domain", attributes={"name":domain})
                        for child in domain_node:
                            logger.info("        %s: %s" %(child.tag, child.text))

                # write out mapping files
                grids = [ ("atm_grid", component_grids[0]), ("lnd_grid", component_grids[1]), ("ocn_grid", component_grids[2]), \
                          ("rof_grid", component_grids[3]), ("glc_grid", component_grids[5]), ("wav_grid", component_grids[6]) ]

                for idx, grid in enumerate(grids):
                    for other_grid in grids[idx+1:]:
                        nodes = self.get_nodes(nodename="gridmap", attributes={grid[0]:grid[1], other_grid[0]:other_grid[1]})
                        for gridmap_node in nodes:
                            for child in gridmap_node:
                                logger.info("    mapping file %s: %s" %(child.tag, child.text))

                logger.info("   ")

                # write out XROT_FLOOD_MODE element TODO: (why is this in there???)
