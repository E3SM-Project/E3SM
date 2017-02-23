#!/usr/bin/env python
"""
config_pes_converter.py -- convert (or verify) config_pesfrom CIME2 format
    to CIME5

The location of these files are needed by the script:
    CIME2: cime/machines-acme/config_pes.xml
    CIME5: config/acme/allactive/config_pesall.xml
"""


from standard_script_setup import *
from CIME.utils import run_cmd
from distutils.spawn import find_executable
import xml.etree.ElementTree as ET
import grid_xml_converter
LOGGER = logging.getLogger(__name__)

###############################################################################
def parse_command_line(args):
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    CIME.utils.setup_standard_logging_options(parser)

    # Set command line options
    parser.add_argument("-cime2file", "--cime2file", help="location of config_grid.xml file in CIME2 repository")
    parser.add_argument("-cime5file", "--cime5file", help="location of config_grids.xml file in CIME5 repository")

    args = parser.parse_args(args[1:])

    CIME.utils.handle_standard_logging_options(args)

    if args.cime2file is None or args.cime5file is None:
        parser.print_help()
        exit()

    return args.cime2file, args.cime5file

class PesNode(grid_xml_converter.DataNode):
    def __init__(self,root):
        self.ignore = False
        super(PesNode, self).__init__(root)

    def __str__(self):
        return ET.tostring(self.xmlnode)

    def setattrib(self, node, tag, key=None):
        if key is None:
            key = tag
        if key in self.data:
            node.set(tag, self.data[key])
        else:
            node.set(tag, 'any')

    def keyvalue(self):
        return "%s:%s:%s:%s" % (self.data['gridname'], self.data['machname'],
                                self.data['pesize'], self.data['compset'])


    def to_cime5(self):
        gridnode = ET.Element('grid')
        self.setattrib(gridnode, 'name', 'gridname')
        machnode = ET.SubElement(gridnode, 'mach')
        self.setattrib(machnode, 'name', 'machname')
        pesnode = ET.SubElement(machnode, 'pes')
        self.setattrib(pesnode, 'compset')
        self.setattrib(pesnode, 'pesize')
        commentnode = ET.SubElement(pesnode, 'comment')
        commentnode.text = "none"
        for d in ['ntasks', 'nthrds', 'rootpe']:
            newnode = ET.SubElement(pesnode, d)
            for comp in ['atm', 'lnd', 'rof', 'ice', 'ocn', 'glc', 'wav', 'cpl']:
                tag = d + '_' + comp
                if tag in self.data[d]:
                    ET.SubElement(newnode, tag).text = str(self.data[d][tag])

        return gridnode



    def __eq__(self, other):
        for k in ['gridname', 'machname', 'pesize', 'compset']:
            if k not in self.data and k not in other.data:
                continue
            if k not in self.data or k not in other.data:
                return False
            if self.data[k] != other.data[k]:
                return False
        for d in ['ntasks', 'nthrds', 'rootpe']:
            for k in self.data[d]:
                if k not in self.data[d] and k not in other.data[d]:
                    continue
                if k not in self.data[d] or k not in other.data[d]:
                    return False
                if self.data[d][k] != other.data[d][k]:
                    return False
        return True

class Cime5PesNode(PesNode):
    def set_data(self, xmlnode):
        for d in ['ntasks', 'nthrds', 'rootpe']:
            self.data[d] = {}
        self.xmlnode = xmlnode
        self.data['gridname'] = xmlnode.get('name')
        machnode = xmlnode.find('mach')
        self.data['machname'] = machnode.get('name')
        pesnode = machnode.find('pes')
        self.data['pesize'] = pesnode.get('pesize')
        self.data['compset'] = pesnode.get('compset')
        commentnode = pesnode.find('comment')
        if commentnode is not None:
            self.data['comment'] = commentnode.text
        for tag in ['ntasks', 'nthrds', 'rootpe']:
            node = pesnode.find(tag)
            for child in node.getchildren():
                self.data[tag][child.tag] = child.text.strip()

class Cime2PesNode(PesNode):
    ISDEFAULT = "-999999"
    DEFAULTS = {'ntasks':'16', 'nthrds':'1', 'rootpe':'0'}
    def set_data(self, xmlnode):
        # Set Defaults
        for d in ['ntasks', 'nthrds', 'rootpe']:
            self.data[d] = {}
        for comp in ['atm', 'lnd', 'ice', 'ocn', 'glc', 'rof', 'wav', 'cpl']:
            self.data['ntasks']['ntasks_' + comp] = self.ISDEFAULT
            self.data['nthrds']['nthrds_' + comp] = self.ISDEFAULT
            self.data['rootpe']['rootpe_' + comp] = self.ISDEFAULT

        # Read in node
        self.xmlnode = xmlnode
        for checktag in ['OS', 'TEST']:
            check = xmlnode.get(checktag)
            if check is not None:
                self.ignore = True
                return
        self.data['machname'] = xmlnode.get('MACH', default='any')
        self.data['gridname'] = xmlnode.get('GRID', default='any')
        self.data['pesize'] = xmlnode.get('PECOUNT', default='any')
        self.data['compset'] = xmlnode.get('CCSM_LCOMPSET', default='any')
        for d in ['ntasks', 'nthrds', 'rootpe']:
            for comp in ['atm', 'lnd', 'ice', 'ocn', 'glc', 'rof', 'wav', 'cpl']:
                tag = d + '_' + comp
                node = xmlnode.find(tag.upper())
                if node is not None:
                    val = node.text.strip()
                    if val[0] == '$':
                        resolvetag = val[1:]
                        if resolvetag == "MAX_TASKS_PER_NODE":
                            val = '-1'
                        else:
                            refnode = xmlnode.find(resolvetag)
                            if refnode is None:
                                # use default value
                                val = self.data[resolvetag.lower()[0:6]][resolvetag.lower()]
                            else:
                                val = xmlnode.find(resolvetag).text.strip()

                    self.data[d][tag] = val
        # Set to defaults. CIME2 had unresolved defaults that referred
        # back to the ATM value, so setting just the ATM value would in effect
        # set all values
        for d in ['ntasks', 'nthrds', 'rootpe']:
            atmtag = d + '_atm'
            if self.data[d][atmtag] == self.ISDEFAULT:
                self.data[d][atmtag] = self.DEFAULTS[d]
            for comp in ['lnd', 'rof', 'ice', 'ocn', 'glc', 'wav', 'cpl']:
                tag = d + '_' + comp
                if self.data[d][tag] == self.ISDEFAULT:
                    self.data[d][tag] = self.data[d][atmtag]




class PesTree(grid_xml_converter.DataTree):
    def __init__(self, xmlfilename):
        # original xml file has bad comments
        import re, StringIO
        with open(xmlfilename,'r') as xmlfile:
            t1 = xmlfile.read()
            t2 = re.sub(r'(?<=<!--)([ -]+)',
                        lambda x: x.group(0).replace('-',' '), t1)
            t3 = re.sub(r'([ -]+)(?=-->)',
                        lambda x: x.group(0).replace('-',' '), t2)
            tempxml = StringIO.StringIO(t3)
            super(PesTree, self).__init__(tempxml)
            tempxml.close()
        
        
    def populate(self):
        xmlnodes = self.root.findall('grid')
        nodeclass = Cime5PesNode

        if len(xmlnodes) == 0:
            xmlnodes = self.root.findall('pes')
            nodeclass = Cime2PesNode
        for xmlnode in xmlnodes:
            datanode = nodeclass(self.root)
            datanode.set_data(xmlnode)
            if not datanode.ignore:
                self.nodes.append(datanode)



    def writexml(self, addlist, newfilename):
        root = ET.Element('config_pes')
        for a, b in addlist:
            if b is not None:
                root.append(ET.Element('REPLACE'))
                root.append(b.to_cime5())
                root.append(ET.Element('WITH'))
            if a is not None:
                root.append(a.to_cime5())
        xmllint = find_executable("xmllint")
        if xmllint is not None:
            run_cmd("%s --format --output %s -"%(xmllint, newfilename),
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
                    oklist.append([anode, bnode])
                else:
                    fixlist.append([anode, bnode])
                break

        if anode in afound:
            continue

        addlist.append([anode, None])



    LOGGER.info("Number of ok nodes: %d" % len(oklist))
    LOGGER.info("Number of wrong nodes: %d" % len(fixlist))
    LOGGER.info("Number of missing nodes: %d" % len(addlist))
    for miss in addlist:
        LOGGER.info(miss[0].keyvalue())
    LOGGER.info("Number of duplicate nodes: %d" % len(duplist))
    for dup in duplist:
        LOGGER.info(dup)
    return [oklist, fixlist, addlist]


def pes_compare():
    cime2file, cime5file = parse_command_line(sys.argv)

    cime2pestree = PesTree(cime2file)
    cime5pestree = PesTree(cime5file)

    LOGGER.info("Comparing config_pes files...")
    oklist, fixlist, addlist = diff_tree(cime2pestree, cime5pestree)
    cime5pestree.postprocess(fixlist, addlist, "tempgrid.xml", cime5file,
                             "badgrid.xml")

if __name__ == "__main__":
    pes_compare()


