"""
Interface to the config_files.xml file.  This class inherits from EntryID.py
"""
from standard_module_setup import *

from generic_xml import GenericXML
from CIME.utils import expect, get_cime_root, get_model

class Testlist(GenericXML):
    def __init__(self,infile):
        """
        initialize an object
        """
        GenericXML.__init__(self,infile)

    def get_tests(self, machine=None, category=None, compiler=None):
        tests = {}
        testnodes = self.get_node("test")
        machatts = {}
        if(machine is not None):
            machatts["name"]      = machine
        if(category is not None):
            machatts["category"]  = category
        if(compiler is not None):
            machatts["compiler"]  = compiler

        for tnode in testnodes:
            machnodes = self.get_node("machine",machatts,root=tnode)

            if(machnodes):
                teststr = "%s.%s.%s"%(tnode.attrib["name"],tnode.attrib["grid"],tnode.attrib["compset"])

                for mach in machnodes:
                    thistest = "%s.%s_%s"%(teststr, mach.attrib["name"],mach.attrib["compiler"])
                    print "Found "+thistest
                    optionnodes = self.get_node("option")


