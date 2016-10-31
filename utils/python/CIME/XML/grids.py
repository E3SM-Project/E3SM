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
        self._version = self.get_version()
        self._comp_gridnames = self._get_grid_names()


    def _get_grid_names(self):
        if self._version == "1.0":
            gridnames = ["atm", "lnd", "ocnice", "rof", "mask", "glc", "wav"]
        elif self._version == "2.0":
            nodes = self.get_nodes("grid")
            gridnames = []
            for node in nodes:
                gn = node.get("name")
                if gn not in gridnames:
                    gridnames.append(gn)
        else:
            expect(False,"Did not recognize config_grids.xml file version")

        return gridnames

    def get_grid_info(self, name, compset):
        """
        Find the matching grid node
        """
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

        # determine component_grids dictionary and grid longname
        lname, component_grids = self._read_config_grids(name, compset)

        gridinfo["GRID"] = lname

        # determine domains given component_grids
        domains  = self._get_domains(component_grids)
        gridinfo.update(domains)

        # determine gridmaps given component_grids
        gridmaps = self._get_gridmaps(component_grids, atmnlev, lndnlev)
        gridinfo.update(gridmaps)

        return gridinfo




    def _read_config_grids(self, name, compset):
        if self._version == "1.0":
            return self._read_config_grids_v1(name, compset)
        elif self._version == "2.0":
            return self._read_config_grids_v2(name, compset)

    def _read_config_grids_v1(self, name, compset):
        """
        read config_grids.xml with v1.0 schema
        """
        component_grids = {}
        nodes = self.get_nodes("grid")
        # first search for all grids that have a compset match - if one is found then return
        for node in nodes:
            if "compset" in node.attrib:
                attrib = node.get("compset")
                compset_match = re.search(attrib,compset)
                if compset_match is not None:
                    alias = self.get_element_text("alias", root=node)
                    sname = self.get_element_text("sname", root=node)
                    lname = self.get_element_text("lname", root=node)
                    if alias == name or lname == name or sname == name:
                        logger.debug("Found node compset match: %s and lname: %s" % (attrib, lname))
                        component_grids = self._get_component_grids(lname)
                        return lname, component_grids
        # if no matches were found for a possible compset match, then search for just a grid match with no
        # compset attribute
        for node in nodes:
            if "compset" not in node.attrib:
                sname = self.get_element_text("sname", root=node)
                alias = self.get_element_text("alias", root=node)
                lname = self.get_element_text("lname", root=node)
                if alias == name or lname == name or sname == name:
                    logger.debug("Found node compset match: %s and lname: %s" % (attrib, lname))
                    component_grids = self._get_component_grids(lname)
                    return lname, component_grids
        expect (False,
                "grid '%s'  is not supported, use manage_case to determine supported grids " %name)


    def _read_config_grids_v2(self, name,compset):
        """
        read config_grids.xml with version 2.0 schema
        """
        component_grids = {}
        model_grid = {}
        for comp_gridname in self._comp_gridnames:
            model_grid[comp_gridname] = None

        # (1) set array of component grid defaults that match current compset

        grid_defaults_node = self.get_node("model_grid_defaults")
        for grid_node in grid_defaults_node:
            name_attrib = grid_node.get("name")
            compset_attrib = grid_node.get("compset")
            compset_match = re.search(compset_attrib, compset)
            if compset_match is not None:
                model_grid[name_attrib] = grid_node.text

        # (2)loop over all of the "model grid" nodes and determine is there an alias match with the
        # input grid name -  if there is an alias match determine if the "compset" and "not_compset"
        # regular expression attributes match the match the input compset

        model_gridnodes = self.get_nodes("model_grid")
        model_gridnode = None
        for node in model_gridnodes:
            found = False
            alias = node.get("alias")
            if alias == name:
                compset_attrib = node.get("compset")
                not_compset_attrib = node.get("not_compset")
                if compset_attrib and not_compset_attrib:
                    compset_match = re.search(compset_attrib, compset)
                    not_compset_match = re.search(not_compset_attrib, compset)
                    if compset_match is not None and not_compset_match is not None:
                        found = True
                        model_gridnode = node
                        logger.debug("Found match for %s with compset_match %s and not_compset_match %s"
                                     % (alias, compset_attrib, not_compset_attrib))
                        break
                elif compset_attrib:
                    compset_match = re.search(compset_attrib, compset)
                    if compset_match is not None:
                        found = True
                        model_gridnode = node
                        logger.debug("Found match for %s with compset_match %s"
                                     % (alias, compset_attrib))
                        break
                elif not_compset_attrib:
                    not_compset_match = re.search(not_compset_attrib, compset)
                    if not_compset_match is None:
                        found = True
                        model_gridnode = node
                        logger.debug("Found match for %s with not_compset_match %s"
                                     % (alias, not_compset_attrib))
                        break
                else:
                    found = True
                    model_gridnode = node
                    logger.debug("Found match for %s" %(alias))
                    break

        # if no match is found in config_grids.xml - exit
        expect(found, "ERROR: no alias was found for %s " %name)

        # for the match - find all of the component grid settings
        grid_nodes = self.get_nodes("grid", root=model_gridnode)
        for grid_node in grid_nodes:
            name = grid_node.get("name")
            value = grid_node.text
            if model_grid[name] != "null":
                model_grid[name] = value
        mask_node = self.get_optional_node("mask",root=model_gridnode)
        if mask_node is not None:
            model_grid["mask"] = mask_node.text
        else:
            model_grid["mask"] = model_grid["ocnice"]

        # determine component grids and associated required domains and gridmaps
        # TODO: this should be in XML, not here
        prefix = {"atm":"a%", "lnd":"l%", "ocnice":"oi%", "rof":"r%", "wav":"w%", "glc":"g%", "mask":"m%"}
        lname = ""
        for component_gridname in self._comp_gridnames:
            if lname:
                lname = lname + "_" + prefix[component_gridname]
            else:
                lname = prefix[component_gridname]
            if model_grid[component_gridname] is not None:
                lname += model_grid[component_gridname]
            else:
                lname += 'null'
        component_grids = self._get_component_grids_from_longname(lname)
        return lname, component_grids

    def _get_domains_v2(self, component_grids):
        """ determine domains dictionary for config_grids.xml v2 schema"""
        # use component_grids to create grids dictionary
        # TODO: this should be in XML, not here
        grids = [("atm", "a%"), ("lnd", "l%"), ("ocn", "o%"), \
                 ("ice", "i%"), ("rof", "r%"), ("glc", "g%"), ("wav", "w%")]
        logger.warn(component_grids)
        domains = {}
        if 'm%' in component_grids:
            mask_name = component_grids['m%']
        for grid in grids:
            grid_name = component_grids[grid[1]]
            domain_node = self.get_optional_node(nodename="domain", attributes={"name":grid_name})
            if domain_node is not None:
                comp_name = grid[0].upper()
                domains[comp_name + "_NX"] = int(self.get_element_text("nx", root=domain_node))
                domains[comp_name + "_NY"] = int(self.get_element_text("ny", root=domain_node))
                domains[comp_name + "_GRID"] = grid_name
                file_name = grid_name + "_DOMAIN_FILE"
                path_name = grid_name + "_DOMAIN_PATH"
                file_nodes = self.get_nodes(nodename="file", root=domain_node)
                for file_node in file_nodes:
                    grid_attrib = file_node.get("grid")
                    mask_attrib = file_node.get("mask")
                    domain_name = ""
                    if grid_attrib is not None and mask_attrib is not None:
                        grid_match = re.search(grid_name.lower(), grid_attrib)
                        mask_match = re.search(mask_name, mask_attrib)
                        if grid_match is not None and mask_match is not None:
                            domain_name = file_node.text
                    elif grid_attrib is not None:
                        grid_match = re.search(grid_name.lower(), grid_attrib)
                        if grid_match is not None:
                            domain_name = file_node.text
                    elif mask_attrib is not None:
                        mask_match = re.search(mask_name.lower(), mask_attrib)
                        if mask_match is not None:
                            domain_name = file_node.text
                    if domain_name:
                        domains[file_name] = os.path.basename(domain_name)
                        domains[path_name] = os.path.dirname(domain_name)
        logger.warn("domains %s"%domains)
        return domains

    def _get_component_grids_from_longname(self, name):
        gridRE = re.compile(r"[_]{0,1}[a-z]{1,2}%")
        grids = gridRE.split(name)[1:]
        prefixes = re.findall("[a-z]+%",name)
        component_grids = {}
        i = 0
        while i < len(grids):
            prefix = prefixes[i]
            grid = grids[i]
            component_grids[prefix] = grid
            i += 1
        component_grids["i%"] = component_grids["oi%"]
        component_grids["o%"] = component_grids["oi%"]
        return component_grids

    def _get_component_grids(self, name):
        gridRE = re.compile(r"[_]{0,1}[a-z]{1,2}%")
        component_grids = gridRE.split(name)[1:]
        return component_grids

    def _get_domains(self, component_grids):
        if self._version == "1.0":
            return self._get_domains_v1(component_grids)
        elif self._version == "2.0":
            return self._get_domains_v2(component_grids)

    def _get_domains_v1(self, component_grids):
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
                domains[grid[0]+"_NX"] = int(self.get_element_text("nx", root=root))
                domains[grid[0]+"_NY"] = int(self.get_element_text("ny", root=root))
                domains[grid[0] + "_GRID"] = grid[1]
                if mask_name is not None:
                    file_ = self.get_element_text("file", attributes={mask_name:mask}, root=root)
                    path  = self.get_element_text("path", attributes={mask_name:mask}, root=root)
                    if file_ is not None:
                        domains[file_name] = file_
                    if path is not None:
                        domains[path_name] = path
        return domains

    def _get_gridmaps(self, component_grids, atmnlev=None, lndnlev=None):
        if self._version == "1.0":
            return self._get_gridmaps_v1(component_grids, atmnlev, lndnlev)
        elif self._version == "2.0":
            return self._get_gridmaps_v2(component_grids)

    def _get_gridmaps_v1(self, component_grids, atmnlev=None, lndnlev=None):
        # mapping files
        grids = [("atm_grid", component_grids[0]), \
                 ("lnd_grid", component_grids[1]), \
                 ("ocn_grid", component_grids[2]), \
                 ("rof_grid", component_grids[3]), \
                 ("glc_grid", component_grids[5]), \
                 ("wav_grid", component_grids[6]), \
                 ("mask"    , component_grids[4])]
        if atmnlev is not None:
            grids[0] += "z"+atmnlev
        if lndnlev is not None:
            grids[1] += "z"+lndnlev

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

        return gridmaps

    def _get_gridmaps_v2(self, component_grids):
        """
        set all mapping files for config_grids.xml v2 schema
        """
        grids = [("atm_grid","a%"), ("lnd_grid","l%"), ("ocn_grid","o%"), \
                 ("rof_grid","r%"), ("glc_grid","g%"), ("wav_grid","w%")]
        gridmaps = {}

        # (1) set all possibly required gridmaps to idmap
        required_gridmaps_node = self.get_node(nodename="required_gridmaps")
        required_gridmap_nodes = self.get_nodes(nodename="required_gridmap", root=required_gridmaps_node)
        for node in required_gridmap_nodes:
            gridmaps[node.text] = "idmap"

        # (2) determine values gridmaps for target grid
        for idx, grid in enumerate(grids):
            for other_grid in grids[idx+1:]:
                gridname = grid[0]
                other_gridname = other_grid[0]
                gridvalue = component_grids[grid[1]]
                other_gridvalue = component_grids[other_grid[1]]
                gridmap_nodes = self.get_nodes(nodename="gridmap",
                                               attributes={gridname:gridvalue, other_gridname:other_gridvalue})
                for gridmap_node in gridmap_nodes:
                    map_nodes = self.get_nodes(nodename="map",root=gridmap_node)
                    for map_node in map_nodes:
                        name = map_node.get("name")
                        value = map_node.text
                        if name is not None and value is not None:
                            gridmaps[name] = value
                            logger.debug(" gridmap name,value are %s: %s" %(name,value))

        # (3) check that all necessary maps are not set to idmap
        griddict = dict(grids)
        for node in required_gridmap_nodes:
            grid1_name = node.get("grid1")
            grid2_name = node.get("grid2")
            prefix1 = griddict[grid1_name]
            prefix2 = griddict[grid2_name]
            grid1_value = component_grids[prefix1]
            grid2_value = component_grids[prefix2]
            if grid1_value is not None and grid2_value is not None:
                if grid1_value != grid2_value and grid1_value != 'null' and grid2_value != 'null':
                    map_ = gridmaps[node.text]
                    if map_ == 'idmap':
                        logger.warning("Warning: missing non-idmap %s for %s, %s and %s %s "
                                       %(node.text, grid1_name, grid1_value, grid2_name, grid2_value))

        return gridmaps

    def print_values(self, long_output=None):
        # write out help message
        helptext = self.get_element_text("help")
        logger.info("%s " %helptext)

        # write out grid elements
        grid_nodes = self.get_nodes(nodename="grid")
        for grid_node in grid_nodes:
            lname = self.get_element_text("lname",root=grid_node)
            sname = self.get_element_text("sname",root=grid_node)
            alias = self.get_element_text("alias",root=grid_node)
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
