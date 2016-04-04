"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
from standard_module_setup import *
import socket
from generic_xml import GenericXML
from files import Files
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Pes(GenericXML):
    def __init__(self, target_component, infile=None, files=None):
        """
        initialize an object
        if a filename is provided it will be used,
        otherwise if a files object is provided it will be used
        otherwise create a files object from default values
        """
        self.pes = None
        self.pes_dir = None
        
        if infile is None:
            if files is None:
                files = Files()
                infile = files.get_value("PES_SPEC_FILE", {"component":target_component})
                self.pes_dir = os.path.dirname(infile)

        GenericXML.__init__(self, infile)

    def find_pes_layout(self, grid, compset, machine, pesize_opts=None):
        pe_select = None
        grid_match = None
        mach_match = None
        compset_match = None
        pesize_match = None

        # Find default which will be overwritten by any match below 
        nodes = self.root.findall(".//grid[@name='any']/mach[@name='any']")
        for node in nodes:
            pes_node = self.get_nodes(nodename="pes", attributes={"pesize":"any", "compset":"any"}, root=node)
            if len(pes_node) == 0:
                continue
            elif len(pes_node) == 1:
                pe_select = node
                break

        # Find nodes with grid= not 'any' and mach='any' - this override the default for $grid_match and $mach_match
        nodes = self.root.findall(".//grid")
        for node in nodes:
            mach_node = self.get_node(nodename="mach", root=node)
            grid_attr = node.attrib["name"]
            mach_attr = mach_node.attrib["name"]
            # go to the next node if grid is 'any' and mach is not 'any'
            if grid_attr == "any" and mach_attr != 'any': 
                continue 
            # determine if there is a match between the grid_attr and the grid longname
            match =  re.search(grid_attr, grid)
            if match:
                grid_match = grid_attr
                mach_match = 'any'
                pe_select = node
                print "1 found grid match for grid_attr %s and  grid %s" %(grid_attr, grid)
                break

        # Find nodes with grid = 'any' and mach = not 'any'
        # The first match found will overwrite the above grid_match and mach_match default
        nodes = self.root.findall(".//grid")
        for node in nodes:
            mach_node = self.get_node(nodename="mach", root=node)
            grid_attr = node.attrib["name"]
            mach_attr = mach_node.attrib["name"]
            # go to the next node if grid is not 'any' and mach is 'any'
            if grid_attr != "any" and mach_attr == 'any': 
                continue 
            # determine if there is a match between the mach_attr and the machine name
            match =  re.search(mach_attr, machine)
            if match:
                grid_match = 'any'
                mach_match = mach_attr
                pe_select = node
                print "2 found machine match for mach_attr %s and  machine %s" %(mach_attr, machine)
                break

        # Find nodes with grid = not 'any' and mach = not 'any'
        # The first match found will overwrite the above grid_match and mach_match default
        nodes = self.root.findall(".//grid")
        for node in nodes:
            mach_node = self.get_node(nodename="mach", root=node)
            grid_attr = node.attrib["name"]
            mach_attr = mach_node.attrib["name"]
            # go to the next node if grid is 'any' and mach is 'any'
            if grid_attr == "any" and mach_attr == 'any': 
                continue 
            # determine if there is a match between the mach_attr and the machine name
            matchm =  re.search(mach_attr, machine)
            matchg =  re.search(grid_attr, grid)
            if matchg and matchm:
                grid_match = grid_attr
                mach_match = mach_attr
                pe_select = node
                print "3 found a grid/machine match for mach_attr %s and  machine %s" %(mach_attr, machine)
                print "3 found a grid/machine match for grid_attr %s and  grid %s" %(grid_attr, grid)
                break

        # Now grid_match and mach_match are set
        # Determine the values of compset_match and pe_size_match
        nodes = self.root.findall(".//grid")
        for node in nodes:
            mach_node = self.get_node(nodename="mach", root=node)
            grid_attr = node.attrib["name"]
            mach_attr = mach_node.attrib["name"]

            # go to the next node if grid and mach do not match grid_match and mach_match, respectively
            matchg = grid_attr == grid_match
            matchm = mach_attr == mach_match
            if not matchg or not matchm:
                continue 

                pes_node = self.get_node(nodename="pes", root=mach_node)
                compset_attr = pes_node.attrib["compset"]
                pesize_attr = pes_node.attrib["pesize"]

                if compset_match is None and pesize_match is None:
                    matchc =  re.search(compset_attr, compset)
                    matchp = None
                    if pesize_opts is not None:
                        matchp = re.search(pesize_attr, pesize_opts)
                    if matchc and matchp:
                        compset_match = compset_attr;
                        pesize_match  = pesize_attr;
                        pe_select = node
                        print "4 found a compset/pesize match for mach_attr %s and  pesize_opts %s" %(mach_attr, pesize_opts)
                        print "4 found a compset/pesize match for compset_attr %s and  compset %s" %(compset_attr, compset)
                        break

                if compset_match is None and pesize_match is None:
                    matchc = re.search(compset_attr, compset)
                    matchp = pesize_attr == 'any'
                    if matchc and matchp:
                        compset_match = compset_attr;
                        pesize_match  = pesize_attr;
                        print "5 found a compset/pesize match for mach_attr %s and  pesize_opts %s" %(mach_attr, pesize_opts)
                        pe_select = node
                        print "5 found a compset/pesize match for compset_attr %s and  compset %s" %(compset_attr, compset)
                        break

                if compset_match is None and pesize_match is None:
                    matchc = compset == 'any'
                    matchp = None
                    if pesize_opts is not None:
                        matchp =  re.search(pesize_attr, pesize_opts)
                    if matchc and matchp:
                        compset_match = compset_attr;
                        pesize_match  = pesize_attr;
                        pesize_match  = pesize_attr;
                        pe_select = node
                        print "6 found a compset/pesize match for mach_attr %s and  pesize_opts %s" %(mach_attr, pesize_opts)
                        print "6 found a compset/pesize match for compset_attr %s and  compset %s" %(compset_attr, compset)
                        break

                if compset_match is None and pesize_match is None:
                    matchc = re.search(compset_attr, compset)
                    matchp = None
                    if pesize_opts is not None:
                        matchp =  re.search(pesize_attr, pesize_opts)
                    if matchc and matchp:
                        compset_match = compset_attr;
                        pesize_match  = pesize_attr;
                        pe_select = node
                        print "7 found a compset/pesize match for mach_attr %s and  pesize_opts %s" %(mach_attr, pesize_opts)
                        print "7 found a compset/pesize match for compset_attr %s and  compset %s" %(compset_attr, compset)
                        break

        pes_ntasks = {}
        nodes = pe_select.findall(".//ntasks/*")
        for node in nodes:
            name = node.tag.upper()
            value = node.text
            pes_ntasks[name] = value
            logger.debug("%s %s " %(name, value))

        pes_nthrds = {}
        nodes = pe_select.findall(".//nthrds/*")
        for node in nodes:
            name = node.tag.upper()
            value = node.text
            pes_nthrds[name] = value
            logger.debug("%s %s " %(name, value))

        pes_rootpe = {}
        nodes = pe_select.findall(".//rootpe/*")
        for node in nodes:
            name = node.tag.upper()
            value = node.text
            pes_rootpe[name] = value
            logger.debug("%s %s " %(name, value))

        logger.info("Pes setting: grid          is %s " %grid)
        logger.info("Pes setting: compset       is %s " %compset)
        logger.info("Pes setting: grid match    is %s " %grid_match )
        logger.info("Pes setting: machine match is %s " %mach_match)
        logger.info("Pes setting: compset_match is %s " %compset_match) 
        logger.info("Pes setting: pesize match  is %s " %pesize_match) 

        return pes_ntasks, pes_nthrds, pes_rootpe
