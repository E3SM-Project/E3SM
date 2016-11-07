"""
Common interface to XML files which follow the grids format,
This is not an abstract class - but inherits from the abstact class GenericXML
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.files import Files
from CIME.XML.generic_xml import GenericXML

logger = logging.getLogger(__name__)

class Grids(GenericXML):

    def __init__(self, infile=None, files=None):
        if infile is None:
            if files is None:
                files = Files()
            infile = files.get_value("GRIDS_SPEC_FILE")
        logger.debug(" Grid specification file is %s" % infile)

        GenericXML.__init__(self, infile)

    def get_grid_info(self, name, compset):
        """
        Find the matching grid node
        """
        nodes = self.get_nodes("grid")
        gridinfo = {}
        atmnlev = None
        lndnlev = None
        #mechanism to specify atm levels
        levmatch = re.match(r"([^_]+)z(\d+)(.*)$", name)
        if  levmatch:
            atmnlev = levmatch.group(2)
            name = levmatch.group(1)+levmatch.group(3)
        #mechanism to specify lnd levels
        levmatch = re.match(r"(.*_)([^_]+)z(\d+)(.*)$", name)
        if  levmatch:
            lndnlev = levmatch.group(3)
            name = levmatch.group(1)+levmatch.group(2)+levmatch.group(4)


        # first search for all grids that have a compset match - if one is found then return
        for node in nodes:
            if "compset" in node.attrib:
                attrib = node.get("compset")
                compset_match = re.search(attrib,compset)
                if compset_match is not None:
                    alias = self.get_value("alias", root=node)
                    sname = self.get_value("sname", root=node)
                    lname = self.get_value("lname", root=node)
                    if alias == name or lname == name or sname == name:
                        logger.debug("Found node compset match: %s and lname: %s" % (attrib, lname))
                        component_grids = self._get_component_grids(lname)
                        domains  = self._get_domains(component_grids)
                        gridmaps = self._get_gridmaps(component_grids, atmnlev, lndnlev)
                        gridinfo.update(domains)
                        gridinfo.update(gridmaps)
                        gridinfo["GRID"] = lname
                        return gridinfo

        # if no matches were found for a possible compset match, then search for just a grid match with no
        # compset attribute
        for node in nodes:
            if "compset" not in node.attrib:
                sname = self.get_value("sname", root=node)
                alias = self.get_value("alias", root=node)
                lname = self.get_value("lname", root=node)
                if alias == name or lname == name or sname == name:
                    logger.debug("Found node compset match: %s and lname: %s" % (attrib, lname))
                    component_grids = self._get_component_grids(lname)
                    gridinfo.update(self._get_domains(component_grids))
                    gridinfo.update(self._get_gridmaps(component_grids, atmnlev, lndnlev))
                    gridinfo["GRID"] = lname
                    return gridinfo

        expect (False,
                "grid '%s'  is not supported, use manage_case to determine supported grids " %name)

    # TODO - API needs to match inherited version
    def get_value(self, item, attributes=None, root=None): # pylint: disable=arguments-differ
        if root is None:
            root = self.root
        node = self.get_optional_node(item,attributes=attributes, root=root)
        if node is not None:
            return node.text

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
        mask_name = None
        for grid in grids:
            file_name = grid[0] + "_DOMAIN_FILE"
            path_name = grid[0] + "_DOMAIN_PATH"
            if grid[0] == "ATM" or grid[0] == "LND":
                mask_name = "lnd_mask"
            if grid[0] == "ICE" or grid[0] == "OCN":
                mask_name = "ocn_mask"
            root = self.get_optional_node(nodename="domain", attributes={"name":grid[1]})
            if root is not None:
                domains[grid[0]+"_NX"] = int(self.get_value("nx", root=root))
                domains[grid[0]+"_NY"] = int(self.get_value("ny", root=root))
                domains[grid[0] + "_GRID"] = grid[1]
                if mask_name is not None:
                    file_ = self.get_value("file", attributes={mask_name:mask}, root=root)
                    path  = self.get_value("path", attributes={mask_name:mask}, root=root)
                    if file_ is not None:
                        domains[file_name] = file_
                    if path is not None:
                        domains[path_name] = path
        return domains

    def _get_gridmaps(self, component_grids, atmnlev=None, lndnlev=None):
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
                            logger.debug(" %s: %s" %gridmap)

        if atmnlev is not None:
            grids[0] += "z"+atmnlev
        if lndnlev is not None:
            grids[1] += "z"+lndnlev

        return gridmaps

    def print_values(self, long_output=None):
        # write out help message
        helptext = self.get_value("help")
        logger.info("%s " %helptext)

        # write out grid elements
        grid_nodes = self.get_nodes(nodename="grid")
        for grid_node in grid_nodes:
            lname = self.get_value("lname",root=grid_node)
            sname = self.get_value("sname",root=grid_node)
            alias = self.get_value("alias",root=grid_node)
            logger.info("------------------------------------------------")
            logger.info("model grid: %s" %lname)
            logger.info("------------------------------------------------")
            if sname is not None:
                logger.info("   short name: %s" %sname)
            if alias is not None:
                logger.info("   alias: %s" %alias)
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
