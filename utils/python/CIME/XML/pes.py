"""
Interface to the config_pes.xml file.  This class inherits from GenericXML.py
"""
from standard_module_setup import *
from generic_xml import GenericXML
from files import Files
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Pes(GenericXML):
    def __init__(self, infile):
        """
        initialize a files object given input pes specification file
        """
        logger.debug("DEBUG: infile is %s"%infile)
        GenericXML.__init__(self, infile)

    def find_pes_layout(self, grid, compset, machine, pesize_opts='M'):
        pe_select = None
        grid_match = None
        mach_match = None
        compset_match = None
        pesize_match = None
        max_points = -1
        # Get all the nodes
        grid_nodes = self.get_nodes("grid")
        for grid_node in grid_nodes:
            grid_match = grid_node.get("name")
            if grid_match == "any" or  re.search(grid_match,grid):
                mach_nodes = self.get_nodes("mach",root=grid_node)
                for mach_node in mach_nodes:
                    mach_match = mach_node.get("name")
                    if mach_match == "any" or re.search(mach_match, machine):
                        pes_nodes = self.get_nodes("pes", root=mach_node)
                        for pes_node in pes_nodes:
                            pesize_match = pes_node.get("pesize")
                            compset_match = pes_node.get("compset")
                            if (pesize_match == "any" or (pesize_opts is not None and \
                                                        re.search(pesize_match, pesize_opts))) and \
                                                        (compset_match == "any" or \
                                                        re.search(compset_match,compset)):

                               points = int(grid_match!="any")*4+int(mach_match!="any")*3+int(compset_match!="any")*2+int(pesize_match!="any")
                               if points > max_points:
                                   pe_select = pes_node
                                   max_points = points
                               elif points == max_points:
                                   expect(False, "We dont expect to be here" )

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
