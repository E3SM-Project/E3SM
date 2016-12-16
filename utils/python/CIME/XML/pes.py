"""
Interface to the config_pes.xml file.  This class inherits from GenericXML.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Pes(GenericXML):

    def __init__(self, infile):
        """
        initialize a files object given input pes specification file
        """
        logger.debug("DEBUG: infile is %s"%infile)
        GenericXML.__init__(self, infile)

    def find_pes_layout(self, grid, compset, machine, pesize_opts='M', mpilib=None):
        pe_select = None
        grid_match = None
        mach_match = None
        compset_match = None
        pesize_match = None
        grid_choice = None
        mach_choice = None
        compset_choice = None
        pesize_choice = None
        other_settings = {}
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

                                points = int(grid_match!="any")*4+int(mach_match!="any")*3+\
                                    int(compset_match!="any")*2+int(pesize_match!="any")
                                if points > max_points:
                                    pe_select = pes_node
                                    max_points = points
                                    mach_choice = mach_match
                                    grid_choice = grid_match
                                    compset_choice = compset_match
                                    pesize_choice = pesize_match
                                elif points == max_points:
                                    logger.warn("mach_choice %s mach_match %s"%(mach_choice, mach_match))
                                    logger.warn("grid_choice %s grid_match %s"%(grid_choice, grid_match))
                                    logger.warn("compset_choice %s compset_match %s"%(compset_choice, compset_match))
                                    logger.warn("pesize_choice %s pesize_match %s"%(pesize_choice, pesize_match))
                                    logger.warn("points = %d"%points)
                                    expect(False, "We dont expect to be here" )

        pes_ntasks, pes_nthrds, pes_rootpe, other_settings = {}, {}, {}, {}
        for node in pe_select:
            vid = node.tag
            logger.debug("vid is %s"%vid)
            if "ntasks" in vid:
                for child in node:
                    pes_ntasks[child.tag.upper()] = child.text
            elif "nthrds" in vid:
                for child in node:
                    pes_nthrds[child.tag.upper()] = child.text
            elif "rootpe" in vid:
                for child in node:
                    pes_rootpe[child.tag.upper()] = child.text
            # if the value is already upper case its something else we are trying to set
            elif vid == node.tag:
                other_settings[vid] = node.text

        logger.info("Pes setting: grid          is %s " %grid)
        logger.info("Pes setting: compset       is %s " %compset)
        logger.info("Pes setting: grid match    is %s " %grid_choice )
        logger.info("Pes setting: machine match is %s " %mach_choice)
        logger.info("Pes setting: compset_match is %s " %compset_choice)
        logger.info("Pes setting: pesize match  is %s " %pesize_choice)
        logger.info("Pes other settings: %s"%other_settings)
        if mpilib == "mpi-serial":
            for i in iter(pes_ntasks):
                pes_ntasks[i] = 1
                pes_rootpe[i] = 1
        return pes_ntasks, pes_nthrds, pes_rootpe, other_settings
