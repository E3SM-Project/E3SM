"""
Interface to the config_files.xml file.  This class inherits from generic_xml.py
It supports versions 1.0 and 2.0 of the testlist.xml file
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

    def _get_testsv1(self, machine=None, category=None, compiler=None):
        tests = []
        compsetnodes = self.get_node("compset")
        machatts = {}
        if(category is not None):
            machatts["testtype"]  = category
        if(compiler is not None):
            machatts["compiler"]  = compiler
        for cnode in compsetnodes:
            gridnodes = self.get_node("grid",root=cnode)
            for gnode in gridnodes:
                testnamenodes = self.get_node("test",root=gnode)
                for tnode in testnamenodes:
                    machnodes = self.get_node("machine",machatts)
                    for mach in machnodes:
                        thistest = {}
                        if (machine is None or (machine is not None and mach.text == machine)):
                            thistest["compiler"] = mach.attrib["compiler"]
                            thistest["category"] = mach.attrib["testtype"]
                            thistest["machine"] = mach.text
                            thistest["testname"] = tnode.attrib["name"]
                            thistest["grid"] = gnode.attrib["name"]
                            thistest["compset"] = cnode.attrib["name"]
                            if ("testmods" in mach.attrib):
                                thistest["testmods"] = mach.attrib["testmods"]
                            tests.append(thistest)
        return tests

    def _get_testsv2(self, machine=None, category=None, compiler=None):
        tests = []
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
            thistest = {}

            if(machnodes):
                for key in tnode.attrib.keys():
                    if(key == "name"):
                        thistest["testname"] = tnode.attrib[key]
                    else:
                        thistest[key] = tnode.attrib[key]
                for mach in machnodes:
                    for key in mach.attrib.keys():
                        if (key == "name"):
                            thistest["machine"] = mach.attrib[key]
                        else:
                            thistest[key] = mach.attrib[key]
                    optionnodes = self.get_node("option", root=mach)
                    for onode in optionnodes:
                       thistest[onode.attrib['name']]=onode.text
                    tests.append(thistest)
        return tests

    def get_tests(self, machine=None, category=None, compiler=None):
        if (self.version == "1.0"):
            return self._get_testsv1(machine,category,compiler)
        elif (self.version == "2.0"):
            return self._get_testsv2(machine,category,compiler)
        else:
            logging.critical("Did not recognize testlist file version %s for file %s"
                             % (self.version,self.filename))

