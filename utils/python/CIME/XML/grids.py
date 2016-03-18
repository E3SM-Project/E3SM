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

    def print_values(self, verbose=None):
        # write out help message
        help_node = self.get_node("help")
        helptext = help_node.text
        print helptext

        # if verbose mode is on, then also obtain nodes for domains and maps
        if verbose is not None:
            domain_nodes = self.get_nodes(nodename="domain")
            gridmap_nodes = self.get_nodes(nodename="gridmap")
            gridRE = re.compile(r"[_]{0,1}[a-z]{1,2}%")

        # write out grid elements 
        grid_nodes = self.get_nodes(nodename="grid")
        for grid_node in grid_nodes:
            lname   = grid_node.find("lname")
            sname   = grid_node.find("sname")
            support = grid_node.find("support")
            alias   = grid_node.find("alias")
            print "model grid: ",lname.text
            if sname is not None:
                print "   short name: ",sname.text
            if alias is not None:
                print "   alias: ",alias.text
            for attr in grid_node.attrib:
                if  attr == 'compset':
                    print "   compset match: ",grid_node.attrib["compset"]

            if verbose is not None:

                # in verbose mode add domain description 
                domain_names = gridRE.split(lname.text)[1:]
                domain_set = list(set(domain_names))
                for domain in domain_set:
                    print "   domain is ",domain
                    domain_node = self.get_node(nodename="domain", attributes={"name":domain}) 
                    for child in domain_node:
                        print  "     ",child.tag,": ",child.text
                    
                # in verbose mode add grid maps 
                # FIXME (mvertens, 3/2016) - how do I make the following more consise and less hard-wired 
                atm_grid  = domain_names[0]
                lnd_grid  = domain_names[1]
                ocn_grid  = domain_names[2]
                rof_grid  = domain_names[3]
                mask_grid = domain_names[4]
                glc_grid  = domain_names[5]
                wav_grid  = domain_names[6]

                nodes_ao = self.get_nodes(nodename="gridmap", attributes={"atm_grid":atm_grid, "ocn_grid":ocn_grid})
                nodes_al = self.get_nodes(nodename="gridmap", attributes={"atm_grid":atm_grid, "lnd_grid":lnd_grid})
                nodes_aw = self.get_nodes(nodename="gridmap", attributes={"atm_grid":atm_grid, "wav_grid":wav_grid})
                nodes_lr = self.get_nodes(nodename="gridmap", attributes={"lnd_grid":lnd_grid, "rof_grid":rof_grid})
                nodes_lg = self.get_nodes(nodename="gridmap", attributes={"lnd_grid":lnd_grid, "glc_grid":glc_grid})
                nodes_ro = self.get_nodes(nodename="gridmap", attributes={"rof_grid":rof_grid, "ocn_grid":ocn_grid})
                nodes_og = self.get_nodes(nodename="gridmap", attributes={"ocn_grid":ocn_grid, "glc_grid":glc_grid})

                for gridmap_node in nodes_ao:
                    for child in gridmap_node:
                        print  "  ",child.tag,": ",child.text
                for gridmap_node in nodes_al:
                    for child in gridmap_node:
                        print  "  ",child.tag,": ",child.text
                for gridmap_node in nodes_aw:
                    for child in gridmap_node:
                        print  "  ",child.tag,": ",child.text
                for gridmap_node in nodes_lr:
                    for child in gridmap_node:
                        print  "  ",child.tag,": ",child.text
                for gridmap_node in nodes_lg:
                    for child in gridmap_node:
                        print  "  ",child.tag,": ",child.text
                for gridmap_node in nodes_ro:
                    for child in gridmap_node:
                        print  "  ",child.tag,": ",child.text
                for gridmap_node in nodes_og:
                    for child in gridmap_node:
                        print  "  ",child.tag,": ",child.text

                print  ""
                
                # write out XROT_FLOOD_MODE element TODO: (why is this in there???)
