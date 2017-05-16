#!/usr/bin/env python
"""
grid_xml_converter.py -- convert (or verify) grid xml from CIME2 format to CIME4
The location of these files are needed by the script:
    CIME2: cime/scripts/Tools/config_grid.xml
    CIME4: cime_config/acme/config_grids.xml
"""

# make sure cime2, cime roots are defined
# use categories
#  GRID CONFIGURATIONS   grid list   domain    grid maps
#    CIME2: cime/scripts/Tools/config_grid.xml
#    CIME4: cime_config/acme/config_grids.xml
#

from standard_script_setup import *
from CIME.utils import run_cmd_no_fail
from distutils.spawn import find_executable
import xml.etree.ElementTree as ET
import operator

logger = logging.getLogger(__name__)

###############################################################################
def parse_command_line(args):
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    CIME.utils.setup_standard_logging_options(parser)

    # Set command line options
    parser.add_argument("-cime2file", "--cime2file", help="location of config_grid.xml file in CIME2 repository")
    parser.add_argument("-cime4file", "--cime4file", help="location of config_grids.xml file in CIME4 repository")

    args = parser.parse_args(args[1:])

    CIME.utils.handle_standard_logging_options(args)

    if args.cime2file is None or args.cime4file is None:
        parser.print_help()
        exit()

    return args.cime2file, args.cime4file



class DataNode(object):
    """
    non-demoninational dictionary of node data:
    """
    def __init__(self, xmlroot):
        self.xmlroot = xmlroot # in case additional information needed
        self.data = {}
        self.name = None
        self.xmlnode = None

    def keyvalue(self):
        return self.data[self.key]

class GridNode(DataNode):
    key = 'lname'
    def __str__(self):
        return ET.tostring(self.xmlnode)


    def to_cime4(self):
        node = ET.Element('grid')
        if 'compset' in self.data and self.data['compset'] is not None:
            node.set('compset',self.data['compset'])

        for k in ['sname','lname','alias','support']:
            if k in self.data and self.data[k] is not None:
                ET.SubElement(node,k).text=self.data[k]

        return node



    def __eq__(self, other):
        for k in ['lname','sname','compset','alias']:
            if k not in self.data and k not in other.data:
                continue
            if k not in self.data or k not in other.data:
                return False
            if self.data[k] != other.data[k]:
                return False
        return True

class Cime2GridNode(GridNode):
    def set_data(self, xmlnode):
        self.xmlnode = xmlnode
        if xmlnode.text is not None:
            self.data['lname'] = xmlnode.text
        for k in ['lname','sname','alias','compset']:
            tmpval = xmlnode.get(k)
            if tmpval is not None:
                self.data[k] = tmpval.strip()
        tmpval = xmlnode.get('support_level')
        if tmpval is not None:
            self.data['support'] = tmpval.strip()

class Cime4GridNode(GridNode):
    def set_data(self, xmlnode):
        self.xmlnode = xmlnode
        for k in ['sname','lname','support','alias']:
            if xmlnode.find(k) is not None:
                self.data[k] = xmlnode.find(k).text.strip()
        if xmlnode.get('compset') is not None:
            self.data['compset'] = xmlnode.get('compset').strip()

class GridmapNode(DataNode):
    def set_data(self, xmlnode):
        self.keys=[]
        self.data['maps']={}
        self.xmlnode = xmlnode
        for k in ['atm_grid','lnd_grid','ocn_grid','rof_grid','glc_grid','wav_grid','ice_grid']:
            att = xmlnode.get(k)
            if att is not None:
                self.data[k]=att.strip()
                self.keys.append(k)
        self.sort()
        for child in xmlnode.getchildren():
            self.data['maps'][child.tag]=child.text.strip()
    def sort(self):
        newlist = sorted(self.keys, key=operator.itemgetter(0))
        self.keys = newlist
    def to_cime4(self):
        node = ET.Element('gridmap')
        for k in ['atm_grid','lnd_grid','ocn_grid','rof_grid','glc_grid']:
            if k in self.data:
                node.set(k,self.data[k])
        for key,value in self.data['maps'].iteritems():
            ET.SubElement(node,key).text=value
        return node
    def __str__(self):
        return str(self.keyvalue()) + str(self.data)
    def __eq__(self,other):
        if self.keyvalue() != other.keyvalue():
            return False
        if len(self.data['maps']) != len(other.data['maps']):
            return False
        for key,value in self.data['maps'].iteritems():
            if key not in other.data['maps'] or value != other.data['maps'][key]:
                return False
        return True

    def keyvalue(self):
        return "%s:%s:%s:%s" % (self.keys[0], self.data[self.keys[0]],
                               self.keys[1], self.data[self.keys[1]])
class DomainNode(DataNode):
    """
    non-demoninational dictionary of domain node information:
    """
    key = 'name'

    def to_cime4(self):
        node = ET.Element('domain')
        node.set('name',self.data['name'])
        for tag in ['nx','ny','desc','support']:
            if tag in self.data:
                ET.SubElement(node,tag).text=self.data[tag]
        for fop in ['file','path']:
            if fop in self.data:
                for comp, mask, filename in self.data[fop]:
                    attribs = {'%s_mask' % comp:mask}
                    ET.SubElement(node,fop,attribs).text=filename
        return node



    def sort(self):
        for fop in ['file','path']:
            newlist = sorted(self.data[fop],key=operator.itemgetter(0))
            self.data[fop] = newlist

    def __eq__(self, other):
        # Check for different name, nx, or ny values
        for k in ['name','nx','ny']:
            if k not in self.data and k not in other.data:
                continue
            if k not in self.data or k not in other.data:
                return False
            if self.data[k] != other.data[k]:
                return False
        # Compare (sorted) file, path lists for equality
        for fop in ['file','path']:
            if fop not in self.data and fop not in other.data:
                contine
            if fop not in self.data or fop not in other.data:
                return False
            if len(self.data[fop]) != len(other.data[fop]):
                return False

            for i in range(0,len(self.data[fop])):
                for j in range(0,2):
                    if self.data[fop][i][j] != other.data[fop][i][j]:
                        return False

        return True

    def __str__(self):
        return str(self.data)

class Cime2DomainNode(DomainNode):
    """
    Read in a domain node from Cime2 xml format
    """
    def set_data(self, xmlnode):
        self.xmlnode = xmlnode
        self.data['name'] = xmlnode.get('name').strip()
        self.data['file'] = []
        self.data['path'] = []
        for tag in ['nx','ny','desc']:
            child = xmlnode.find(tag)
            if child is not None:
                self.data[tag] = child.text

        # Find any griddom nodes that match this name
        griddoms = self.xmlroot.findall('.griddom[@grid="%s"]' % self.data['name'])
        for gd in griddoms:
            mask = gd.get('mask')
            for comp in ['ATM','LND','OCN','ICE']:
                for fop in ['FILE','PATH']:
                    tag = '%s_DOMAIN_%s' % (comp,fop)
                    n = gd.find(tag)
                    if n is not None:
                        self.data[fop.lower()].append([comp.lower(),mask,n.text])
        # sort the file and path entries
        self.sort()

class Cime4DomainNode(DomainNode):
    """
    Read in a domain node from Cime4 xml format
    """
    def set_data(self, xmlnode):
        self.xmlnode = xmlnode
        self.data['name'] = xmlnode.get('name')
        self.data['file'] = []
        self.data['path'] = []
        for tag in ['nx','ny','desc','support']:
            child = xmlnode.find(tag)
            if child is not None:
                self.data[tag] = child.text
        for comp in ['lnd','atm','ocn','ice']:
            masktag = '%s_mask' % comp
            for fop in ['file','path']:
                fopnodes = xmlnode.findall('%s[@%s]' % (fop,masktag))
                for n in fopnodes:
                    mask = n.get(masktag)
                    filename = n.text.strip()
                    self.data[fop].append([comp,mask,filename])

        # sort the file and path entries
        self.sort()


class DataTree(object):
    def __init__(self, xmlfilename):
        self.xmlfilename = xmlfilename
        self.doc = ET.parse(xmlfilename)
        self.root = self.doc.getroot()
        self.index = 0
        self.n = 0
        self.nodes = []
        self.populate()

    def populate(self,expression,classname):
        xmlnodes = self.root.findall(expression)
        for xmlnode in xmlnodes:
            datanode = classname(self.root)
            datanode.set_data(xmlnode)
            self.nodes.append(datanode)

    def next(self):
        if self.index >= len(self.nodes):
            self.index=0
            raise StopIteration
        if self.index < len(self.nodes):
            self.index += 1
            return self.nodes[self.index-1]

    def __iter__(self):
        return self

    def postprocess(self, fixlist, addlist, newxmlfile, badxmlfile):
        if len(addlist) > 0:
            logger.info("\n\nWriting suggested nodes to %s" % newxmlfile)
            logger.info("Copy 'grid' nodes into corresponding location in")
            logger.info("cime_config/acme/config_grids.xml")
            self.writexml(addlist,newxmlfile)
            self.writexml(fixlist,badxmlfile)
            if len(fixlist) > 0:
                logger.info("Some nodes should be removed from")
                logger.info("cime_config/acme/config_grids.xml. These nodes")
                logger.info("have been written to %s" % badxmlfile)

class GridTree(DataTree):
    def populate(self):
        xmlnodes = self.root.findall('GRID')
        nodeclass=Cime2GridNode
        if len(xmlnodes) == 0:
            xmlnodes = self.root.findall('./grids/grid')
            nodeclass=Cime4GridNode

        for xmlnode in xmlnodes:
            datanode = nodeclass(self.root)
            datanode.set_data(xmlnode)
            self.nodes.append(datanode)

    def writexml(self,addlist,newfilename):
        root = ET.Element('grid_data')
        grids = ET.SubElement(root,'grids')
        for a, b in addlist:
            if b is not None:
                grids.append(ET.Element('REPLACE'))
                grids.append(b.to_cime4())
                grids.append(ET.Element('WITH'))

            if a is not None:
                grids.append(a.to_cime4())
        xmllint = find_executable("xmllint")
        if xmllint is not None:
            run_cmd_no_fail("%s --format --output %s -"%(xmllint,newfilename),
                    input_str=ET.tostring(root))


class DomainTree(DataTree):
    def populate(self):
        xmlnodes = self.root.findall('gridhorz')
        nodeclass=Cime2DomainNode
        if len(xmlnodes) == 0:
            xmlnodes = self.root.findall('./domains/domain')
            nodeclass=Cime4DomainNode

        for node in xmlnodes:
            datanode = nodeclass(self.root)
            datanode.set_data(node)
            self.nodes.append(datanode)

    def writexml(self,addlist,newfilename):
        root = ET.Element('grid_data')
        domains = ET.SubElement(root,'domains')
        for a, b in addlist:
            if b is not None:
                domains.append(ET.Element('REPLACE'))
                domains.append(b.to_cime4())
                domains.append(ET.Element('WITH'))
            if a is not None:
                domains.append(a.to_cime4())
        xmllint = find_executable("xmllint")
        if xmllint is not None:
            run_cmd_no_fail("%s --format --output %s -"%(xmllint,newfilename),
                    input_str=ET.tostring(root))

class GridmapTree(DataTree):
    def populate(self):
        xmlnodes = self.root.findall('gridmap')
        if len(xmlnodes) == 0:
            xmlnodes = self.root.findall('./gridmaps/gridmap')
        for xmlnode in xmlnodes:
            datanode = GridmapNode(self.root)
            datanode.set_data(xmlnode)
            self.nodes.append(datanode)

    def writexml(self,addlist,newfilename):
        root = ET.Element('gridmaps')
        gridmaps = ET.SubElement(root,'gridmap')
        for a, b in addlist:
            if b is not None:
                gridmaps.append(ET.Element('REPLACE'))
                gridmaps.append(b.to_cime4())
                gridmaps.append(ET.Element('WITH'))
            if a is not None:
                gridmaps.append(a.to_cime4())
        xmllint = find_executable("xmllint")
        if xmllint is not None:
            run_cmd_no_fail("%s --format --output %s -"%(xmllint,newfilename),
                    input_str=ET.tostring(root))

def diff_tree(atree, btree):
    afound = []
    bfound = []
    oklist = []
    fixlist = []
    addlist = []
    duplist = []
    bkeys = []
    for bnode in btree.nodes:
        if bnode.keyvalue() in bkeys:
            duplist.append(bnode.keyvalue())
        else:
            bkeys.append(bnode.keyvalue())


    for anode in atree.nodes:
        for bnode in btree.nodes:
            if bnode in bfound:
                continue
            if anode.keyvalue() == bnode.keyvalue():
                afound.append(anode)
                bfound.append(bnode)

                if anode == bnode:
                    oklist.append([anode,bnode])
                else:
                    fixlist.append([anode,bnode])
                break

        if anode in afound:
            continue

        addlist.append([anode,None])



    logger.info("Number of ok nodes: %d" % len(oklist))
    logger.info("Number of wrong nodes: %d" % len(fixlist))
    logger.info("Number of missing nodes: %d" % len(addlist))
    logger.info("Number of duplicate nodes: %d" % len(duplist))
    for dup in duplist:
        logger.info(dup)
    return [oklist,fixlist,addlist]




def grid_compare():
    cime2file, cime4file = parse_command_line(sys.argv)

    cime2gridtree = GridTree(cime2file)
    cime4gridtree = GridTree(cime4file)
    cime2domaintree = DomainTree(cime2file)
    cime4domaintree = DomainTree(cime4file)
    cime2gridmaptree = GridmapTree(cime2file)
    cime4gridmaptree = GridmapTree(cime4file)

    logger.info("Comparing grid nodes...")
    oklist, fixlist, addlist = diff_tree(cime2gridtree, cime4gridtree)
    cime4gridtree.postprocess(fixlist, addlist, "tempgrid.xml","badgrid.xml")

    oklist, fixlist, addlist = diff_tree(cime2domaintree, cime4domaintree)
    cime4domaintree.postprocess(fixlist, addlist, "tempdomain.xml","baddomain.xml")

    oklist, fixlist, addlist = diff_tree(cime2gridmaptree, cime4gridmaptree)
    cime4gridmaptree.postprocess(fixlist, addlist, "tempgridmap.xml","badgridmap.xml")

if __name__=="__main__":
    grid_compare()
