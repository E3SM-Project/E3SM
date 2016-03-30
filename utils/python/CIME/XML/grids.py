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

    def get_grid_info(self, name, compset):
        # FIXME - the following can be simplified
        # search for all grid nodes
        nodes = self.get_nodes("grid")

        # first search for all grids that have a compset match - if one is found then return
        for node in nodes:
            if node.attrib:
                attrib = node.attrib["compset"]
                compsetRE = re.compile(attrib)
                compset_match = re.search(attrib,compset)
                if compset_match is not None:
                    alias = self.get_node("alias",root=node)
                    sname = self.get_node("sname",root=node)
                    lname = self.get_node("lname",root=node)
                    if alias.text == name or lname.text == name or sname.text == name:
                        logger.debug("Found node compset match: %s and lname: %s" % (attrib, lname.text))
                        component_grids = self._get_component_grids(lname.text)
                        domains  = self._get_domains(component_grids)
                        gridmaps = self._get_gridmaps(component_grids)
                        gridinfo = {}
                        gridinfo.update(domains)
                        gridinfo.update(gridmaps)
                        gridinfo["GRID"] = lname.text
                        return gridinfo


        # if no matches were found for a possible compset match, then search for just a grid match with no
        # compset attribute
        for node in nodes:
            if not node.attrib:
                sname = self.get_node("sname",root=node)
                alias = self.get_node("alias",root=node)
                lname = self.get_node("lname",root=node)
                if alias.text == name or lname.text == name or sname.text == name:
                    logger.debug("Found node compset match: %s and lname: %s" % (attrib, lname.text))
                    component_grids = self._get_component_grids(lname.text)
                    domain_file, domain_path = self._get_domains(component_grids)
                    gridmaps = self._get_gridmaps(component_grids)
                    return longname.text, component_grids, domain_file, domain_path, gridmaps

    def _get_component_grids(self, name):
        gridRE = re.compile(r"[_]{0,1}[a-z]{1,2}%")
        component_grids = gridRE.split(name)[1:]

        return component_grids

    def _get_domains(self, component_grids):
        # use component_grids to create grids dictionary
        grids = [("ATM", component_grids[0]), \
                 ("LND", component_grids[1]), \
                 ("OCN", component_grids[2]), \
                 ("ICE", component_grids[2]), \
                 ("ROF", component_grids[3]), \
                 ("GLC", component_grids[5]), \
                 ("WAV", component_grids[6])]
        mask = component_grids[4]

        domains = {}
        for idx, grid in enumerate(grids):
            file_name = grid[0] + "_DOMAIN_FILE"
            path_name = grid[0] + "_DOMAIN_PATH"
            mask_name = None
            if grid[0] == "ATM" or grid[0] == "LND":
                mask_name = "lnd_mask"
            if grid[0] == "ICE" or grid[0] == "OCN":
                mask_name = "ocn_mask"
            root = self.get_optional_node(nodename="domain", attributes={"name":grid[1]})
            if root is not None:
                domains[grid[0]+"_NX"] = int(self.get_optional_node(nodename="nx", root=root).text)
                domains[grid[0]+"_NY"] = int(self.get_optional_node(nodename="ny", root=root).text)
                domains[grid[0] + "_GRID"] = grid[1]
                if mask_name is not None:
                    file = self.get_optional_node(nodename="file", attributes={mask_name:mask}, root=root)
                    path = self.get_optional_node(nodename="path", attributes={mask_name:mask}, root=root)
                    if file is not None:
                        domains[file_name] = file.text
                        if path is not None:
                            domains[path_name] = path.text

        return domains

    def _get_gridmaps(self, component_grids):
        # mapping files
        grids = [("atm_grid", component_grids[0]), \
                 ("lnd_grid", component_grids[1]), \
                 ("ocn_grid", component_grids[2]), \
                 ("rof_grid", component_grids[3]), \
                 ("glc_grid", component_grids[5]), \
                 ("wav_grid", component_grids[6]), \
                 ("mask"    , component_grids[4])]

        gridmaps = {}
        for idx, grid in enumerate(grids):
            for other_grid in grids[idx+1:]:
                nodes = self.get_nodes(nodename="gridmap", attributes={grid[0]:grid[1], other_grid[0]:other_grid[1]})
                for gridmap_node in nodes:
                    for child in gridmap_node:
                        gridmap = (child.tag, child.text)
                        if gridmap is not None:
                            gridmaps[child.tag] = child.text
                            logger.info(" %s: %s" %gridmap)

        return gridmaps


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

