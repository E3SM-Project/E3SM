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
            testnodes = self.get_children('test')
            for node in testnodes:
                self._testnodes[self.get(node, "name")] = node

    def set_header(self, testroot, machine, testid, baselinetag=None, baselineroot=None):
        tlelem = self.make_child("testlist")

        for name, text in [ ("testroot", testroot), ("machine", machine), ("testid", testid), ("baselinetag", baselinetag), ("baselineroot", baselineroot) ]:
            if text is not None:
                self.make_child(name, root=tlelem, text=text)

        self._testlist_node = tlelem

    def add_test(self, compiler, mpilib, testname):
        expect(testname not in self._testnodes, "Test {} already in testlist".format(testname))

        telem = self.make_child("test", attributes={"name":testname}, root=self._testlist_node)

        for name, text in [ ("compiler", compiler), ("mpilib", mpilib) ]:
            self.make_child(name, root=telem, text=text)

        self._testnodes[testname] = telem

    def update_test_status(self, testname, phase, status):
        expect(testname in self._testnodes, "Test {} not defined in testlist".format(testname))
        root = self._testnodes[testname]
        pnode = self.get_optional_child("section", {"name":phase}, root=root)
        if pnode is not None:
            self.set(pnode, "status", status)
        else:
            self.make_child("section", {"name":phase, "status":status}, root=root)
