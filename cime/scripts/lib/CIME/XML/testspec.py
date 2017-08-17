"""
Interface to the testspec.xml file.  This class inherits from generic_xml.py
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.generic_xml import GenericXML

logger = logging.getLogger(__name__)

class TestSpec(GenericXML):

    def __init__(self, infile):
        """
        initialize an object
        """
        GenericXML.__init__(self, infile)
        self._testnodes = {}
        self._testlist_node = None
        if os.path.isfile(infile):
            testnodes = self.get_nodes('test')
            for node in testnodes:
                self._testnodes[node.get("name")] = node

    def set_header(self, testroot, machine, testid, baselinetag=None, baselineroot=None):
        tlelem = ET.Element('testlist')
        elem = ET.Element('testroot')
        elem.text = testroot
        tlelem.append(elem)
        elem = ET.Element('machine')
        elem.text = machine
        tlelem.append(elem)
        elem = ET.Element('testid')
        elem.text = testid
        tlelem.append(elem)
        if baselinetag is not None:
            elem = ET.Element('baselinetag')
            elem.text = baselinetag
            tlelem.append(elem)
        if baselineroot is not None:
            elem = ET.Element('baselineroot')
            elem.text = baselineroot
            tlelem.append(elem)
        self.add_child(tlelem)
        self._testlist_node = tlelem

    def add_test(self, compiler, mpilib, testname):
        expect(testname not in self._testnodes, "Test {} already in testlist".format(testname))
        telem = ET.Element("test")
        telem.set("name", testname)
        elem = ET.Element("compiler")
        elem.text = compiler
        telem.append(elem)
        elem = ET.Element("mpilib")
        elem.text = mpilib
        telem.append(elem)
        self.add_child(telem, root=self._testlist_node)
        self._testnodes[testname] = telem

    def update_test_status(self, testname, phase, status):
        expect(testname in self._testnodes, "Test {} not defined in testlist".format(testname))
        root = self._testnodes[testname]
        pnode = self.get_optional_node("section", {"name":phase}, root=root)
        if pnode is not None:
            pnode.set("status", status)
        else:
            pnode = ET.Element("section")
            pnode.set("name", phase)
            pnode.set("status", status)
            root.append(pnode)




