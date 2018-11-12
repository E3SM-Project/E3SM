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
        if files is None:
            files = Files()
        if infile is None:
            infile = files.get_value("GRIDS_SPEC_FILE")
        logger.debug(" Grid specification file is {}".format(infile))
        schema = files.get_schema("GRIDS_SPEC_FILE")

        GenericXML.__init__(self, infile, schema)
        self._version = self.get_version()

        self._comp_gridnames = self._get_grid_names()

    def _get_grid_names(self):
        grids = self.get_child("grids")
        model_grid_defaults = self.get_child("model_grid_defaults", root=grids)
        nodes = self.get_children("grid", root=model_grid_defaults)
        gridnames = []
        for node in nodes:
            gn = self.get(node, "name")
            if gn not in gridnames:
                gridnames.append(gn)
        if "mask" not in gridnames:
            gridnames.append("mask")

        return gridnames

    def get_grid_info(self, name, compset, driver):
        """
        Find the matching grid node
        """
        gridinfo = {}
        atmnlev = None
        lndnlev = None

        #mechanism to specify atm levels
        atmlevregex = re.compile(r"([^_]+)z(\d+)(.*)$")
        levmatch = re.match(atmlevregex, name)
        if  levmatch:
            atmnlev = levmatch.group(2)
            name = levmatch.group(1)+levmatch.group(3)

        #mechanism to specify lnd levels
        lndlevregex = re.compile(r"(.*_)([^_]+)z(\d+)(_[^m].*)$")
        levmatch = re.match(lndlevregex, name)
        if  levmatch:
            lndnlev = levmatch.group(3)
            name = levmatch.group(1)+levmatch.group(2)+levmatch.group(4)

        # determine component_grids dictionary and grid longname
        lname, component_grids = self._read_config_grids(name, compset, atmnlev, lndnlev)
        gridinfo["GRID"] = lname

        # determine domains given component_grids
        domains  = self._get_domains(component_grids, atmlevregex, lndlevregex, driver)

        gridinfo.update(domains)

        # determine gridmaps given component_grids
        gridmaps = self._get_gridmaps(component_grids)
        gridinfo.update(gridmaps)

        return gridinfo

    def _read_config_grids(self, name, compset, atmnlev, lndnlev):
        """
        read config_grids.xml with version 2.0 schema
        """
        component_grids = {}
        model_grid = {}
        for comp_gridname in self._comp_gridnames:
            model_grid[comp_gridname] = None

        # (1) set array of component grid defaults that match current compset
        grids_node = self.get_child("grids")
        grid_defaults_node = self.get_child("model_grid_defaults", root=grids_node)
        for grid_node in self.get_children("grid", root=grid_defaults_node):
            name_attrib = self.get(grid_node, "name")
            compset_attrib = self.get(grid_node, "compset")
            compset_match = re.search(compset_attrib, compset)
            if compset_match is not None:
                model_grid[name_attrib] = self.text(grid_node)

        # (2)loop over all of the "model grid" nodes and determine is there an alias match with the
        # input grid name -  if there is an alias match determine if the "compset" and "not_compset"
        # regular expression attributes match the match the input compset

        model_gridnodes = self.get_children("model_grid", root=grids_node)
        model_gridnode = None
        foundalias = False
        for node in model_gridnodes:
            alias = self.get(node, "alias")
            if alias == name:
                foundalias = True
                foundcompset = False
                compset_attrib = self.get(node, "compset")
                not_compset_attrib = self.get(node, "not_compset")
                if compset_attrib and not_compset_attrib:
                    compset_match = re.search(compset_attrib, compset)
                    not_compset_match = re.search(not_compset_attrib, compset)
                    if compset_match is not None and not_compset_match is not None:
                        foundcompset = True
                        model_gridnode = node
                        logger.debug("Found match for {} with compset_match {} and not_compset_match {}"
                                     .format(alias, compset_attrib, not_compset_attrib))
                        break
                elif compset_attrib:
                    compset_match = re.search(compset_attrib, compset)
                    if compset_match is not None:
                        foundcompset = True
                        model_gridnode = node
                        logger.debug("Found match for {} with compset_match {}"
                                     .format(alias, compset_attrib))
                        break
                elif not_compset_attrib:
                    not_compset_match = re.search(not_compset_attrib, compset)
                    if not_compset_match is None:
                        foundcompset = True
                        model_gridnode = node
                        logger.debug("Found match for {} with not_compset_match {}"
                                     .format(alias, not_compset_attrib))
                        break
                else:
                    foundcompset = True
                    model_gridnode = node
                    logger.debug("Found match for {}".format(alias))
                    break
        expect(foundalias, "no alias {} defined".format(name))
        # if no match is found in config_grids.xml - exit
        expect(foundcompset, "grid alias {} not valid for compset {}".format(name, compset))

        # for the match - find all of the component grid settings
        grid_nodes = self.get_children("grid", root=model_gridnode)
        for grid_node in grid_nodes:
            name = self.get(grid_node, "name")
            value = self.text(grid_node)
            if model_grid[name] != "null":
                model_grid[name] = value
        mask_node = self.get_optional_child("mask",root=model_gridnode)
        if mask_node is not None:
            model_grid["mask"] = self.text(mask_node)
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
                if component_gridname == 'atm' and atmnlev is not None:
                    if not ("a{:n}ull" in lname):
                        lname += "z" + atmnlev

                elif component_gridname == 'lnd' and lndnlev is not None:
                    if not ("l{:n}ull" in lname):
                        lname += "z" + lndnlev

            else:
                lname += 'null'
        component_grids = self._get_component_grids_from_longname(lname)
        return lname, component_grids

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

    def _get_domains(self, component_grids, atmlevregex, lndlevregex, driver):
        """ determine domains dictionary for config_grids.xml v2 schema"""
        # use component_grids to create grids dictionary
        # TODO: this should be in XML, not here
        grids = [("atm", "a%"), ("lnd", "l%"), ("ocn", "o%"), ("mask", "m%"),\
                 ("ice", "i%"), ("rof", "r%"), ("glc", "g%"), ("wav", "w%")]
        domains = {}
        mask_name = None
        if 'm%' in component_grids:
            mask_name = component_grids['m%']
        else:
            mask_name = component_grids['oi%']

        for grid in grids:
            grid_name = component_grids[grid[1]]

            # Determine grid name with no nlev suffix if there is one
            grid_name_nonlev = grid_name
            levmatch = re.match(atmlevregex, grid_name)
            if  levmatch:
                grid_name_nonlev = levmatch.group(1)+levmatch.group(3)
            levmatch = re.match(lndlevregex, grid_name)
            if  levmatch:
                grid_name_nonlev = levmatch.group(1)+levmatch.group(2)+levmatch.group(4)

            # Determine all domain information search for the grid name with no level suffix in config_grids.xml
            domain_node = self.get_optional_child("domain", attributes={"name":grid_name_nonlev},
                                                  root=self.get_child("domains"))
            if domain_node is not None:
                comp_name = grid[0].upper()

                # determine xml variable name
                if not comp_name == "MASK":
                    domains[comp_name + "_NX"] = int(self.get_element_text("nx", root=domain_node))
                    domains[comp_name + "_NY"] = int(self.get_element_text("ny", root=domain_node))
                    file_name = comp_name + "_DOMAIN_FILE"
                    path_name = comp_name + "_DOMAIN_PATH"
                    mesh_name = comp_name + "_DOMAIN_MESH"

                # set up dictionary of domain files for every component
                domains[comp_name + "_GRID"] = grid_name

                file_nodes = self.get_children("file", root=domain_node)
                for file_node in file_nodes:
                    grid_attrib = self.get(file_node, "grid")
                    mask_attrib = self.get(file_node, "mask")
                    domain_name = ""
                    if grid_attrib is not None and mask_attrib is not None:
                        grid_match = re.search(comp_name.lower(), grid_attrib)
                        mask_match = False
                        if mask_name is not None:
                            mask_match = mask_name == mask_attrib
                        if grid_match is not None and mask_match:
                            domain_name = self.text(file_node)
                    elif grid_attrib is not None:
                        grid_match = re.search(comp_name.lower(), grid_attrib)
                        if grid_match is not None:
                            domain_name = self.text(file_node)
                    elif mask_attrib is not None:
                        mask_match = mask_name == mask_attrib
                        if mask_match:
                            domain_name = self.text(file_node)
                    if domain_name:
                        domains[file_name] = os.path.basename(domain_name)
                        path = os.path.dirname(domain_name)
                        if len(path) > 0:
                            domains[path_name] = path

                if not comp_name == "MASK":
                    mesh_nodes = self.get_children("mesh", root=domain_node)
                    for mesh_node in mesh_nodes:
                        driver_attrib = self.get(mesh_node, "driver")
                        if driver == driver_attrib:
                            domains[mesh_name] = self.text(mesh_node)

        return domains

    def _get_gridmaps(self, component_grids):
        """
        set all mapping files for config_grids.xml v2 schema
        """
        grids = [("atm_grid","a%"), ("lnd_grid","l%"), ("ocn_grid","o%"), \
                 ("rof_grid","r%"), ("glc_grid","g%"), ("wav_grid","w%")]
        gridmaps = {}

        # (1) set all possibly required gridmaps to idmap
        required_gridmaps_node = self.get_child("required_gridmaps")
        required_gridmap_nodes = self.get_children("required_gridmap", root=required_gridmaps_node)
        for node in required_gridmap_nodes:
            gridmaps[self.text(node)] = "idmap"

        # (2) determine values gridmaps for target grid
        for idx, grid in enumerate(grids):
            for other_grid in grids[idx+1:]:
                gridname = grid[0]
                other_gridname = other_grid[0]
                gridvalue = component_grids[grid[1]]
                if gridname == "atm_grid":
                    atm_gridvalue = gridvalue
                other_gridvalue = component_grids[other_grid[1]]
                gridmap_nodes = self.get_children("gridmap", root=self.get_child("gridmaps"),
                                               attributes={gridname:gridvalue, other_gridname:other_gridvalue})
                for gridmap_node in gridmap_nodes:
                    expect(len(self.attrib(gridmap_node)) == 2,
                           " Bad attribute count in gridmap node %s"%self.attrib(gridmap_node))
                    map_nodes = self.get_children("map",root=gridmap_node)
                    for map_node in map_nodes:
                        name = self.get(map_node, "name")
                        value = self.text(map_node)
                        if name is not None and value is not None:
                            gridmaps[name] = value
                            logger.debug(" gridmap name,value are {}: {}"
                                         .format(name,value))

        # (3) check that all necessary maps are not set to idmap
        griddict = dict(grids)
        for node in required_gridmap_nodes:
            grid1_name = self.get(node, "grid1")
            grid2_name = self.get(node, "grid2")
            prefix1 = griddict[grid1_name]
            prefix2 = griddict[grid2_name]
            grid1_value = component_grids[prefix1]
            grid2_value = component_grids[prefix2]
            if grid1_value is not None and grid2_value is not None:
                if grid1_value != grid2_value and grid1_value != 'null' and grid2_value != 'null':
                    map_ = gridmaps[self.text(node)]
                    if map_ == 'idmap':
                        if grid1_name == "ocn_grid" and grid1_value == atm_gridvalue:
                            logger.debug('ocn_grid == atm_grid so this is not an idmap error')
                        else:
                            logger.warning("Warning: missing non-idmap {} for {}, {} and {} {} ".format(self.text(node), grid1_name, grid1_value, grid2_name, grid2_value))

        return gridmaps

    def print_values(self, long_output=None):
        # write out help message
        helptext = self.get_element_text("help")
        logger.info("{} ".format(helptext))

        logger.info("{:5s}-------------------------------------------------------------".format(""))
        logger.info("{:10s}  default component grids:\n".format(""))
        logger.info("     component         compset       value " )
        logger.info("{:5s}-------------------------------------------------------------".format(""))
        default_nodes = self.get_children("model_grid_defaults", root=self.get_child("grids"))
        for default_node in default_nodes:
            grid_nodes = self.get_children("grid", root=default_node)
            for grid_node in grid_nodes:
                name = self.get(grid_node, "name")
                compset = self.get(grid_node, "compset")
                value = self.text(grid_node)
                logger.info("     {:6s}   {:15s}   {:10s}".format(name, compset, value))
        logger.info("{:5s}-------------------------------------------------------------".format(""))

        domains = {}
        if long_output is not None:
            domain_nodes = self.get_children("domain",root=self.get_child("domains"))
            for domain_node in domain_nodes:
                name = self.get(domain_node, 'name')
                if name == 'null':
                    continue
                desc = self.text(self.get_child("desc", root=domain_node))
                files = ""
                file_nodes = self.get_children("file", root=domain_node)
                for file_node in file_nodes:
                    filename = self.text(file_node)
                    mask_attrib = self.get(file_node, "mask")
                    grid_attrib = self.get(file_node, "grid")
                    files += "\n       " + filename
                    if mask_attrib or grid_attrib:
                        files += " (only for"
                    if mask_attrib:
                        files += " mask: " + mask_attrib
                    if grid_attrib:
                        files += " grid match: " + grid_attrib
                    if mask_attrib or grid_attrib:
                        files += ")"
                domains[name] = "\n       {} with domain file(s): {} ".format(desc, files)

        model_grid_nodes = self.get_children("model_grid", root=self.get_child("grids"))
        for model_grid_node in model_grid_nodes:
            alias = self.get(model_grid_node, "alias")
            compset = self.get(model_grid_node, "compset")
            not_compset = self.get(model_grid_node, "not_compset")
            restriction = ""
            if compset:
                restriction += "only for compsets that are {} ".format(compset)
            if not_compset:
                restriction += "only for compsets that are not {} ".format(not_compset)
            if restriction:
                logger.info("\n     alias: {} ({})".format(alias,restriction))
            else:
                logger.info("\n     alias: {}".format(alias))
            grid_nodes = self.get_children("grid", root=model_grid_node)
            grids = ""
            gridnames = []
            for grid_node in grid_nodes:
                gridnames.append(self.text(grid_node))
                grids += self.get(grid_node, "name") + ":" + self.text(grid_node) + "  "
            logger.info("       non-default grids are: {}".format(grids))
            mask_nodes = self.get_children("mask", root=model_grid_node)
            for mask_node in mask_nodes:
                logger.info("       mask is: {}".format(self.text(mask_node)))
            if long_output is not None:
                gridnames = set(gridnames)
                for gridname in gridnames:
                    if gridname != "null":
                        logger.info ("    {}".format(domains[gridname]))
