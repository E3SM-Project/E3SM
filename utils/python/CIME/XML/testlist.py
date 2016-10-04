"""
Interface to the config_files.xml file.  This class inherits from generic_xml.py
It supports versions 1.0 and 2.0 of the testlist.xml file

In version 2 of the file options can be specified to further refine a test or set of tests
they can be specified at the <machines> level or at the level of a particular machine.
currently supported options are walltime which sets the wallclock time for a given test
in the queing system, and memleak_tolerance which specifies the relative memory growth expected for a given test.

"""
from CIME.XML.standard_module_setup import *

from CIME.XML.generic_xml import GenericXML

logger = logging.getLogger(__name__)

class Testlist(GenericXML):

    def __init__(self,infile):
        """
        initialize an object
        """
        GenericXML.__init__(self,infile)

    def _get_testsv1(self, machine=None, category=None, compiler=None):
        tests = []

        machatts = {}
        if category is not None:
            machatts["testtype"] = category
        if compiler is not None:
            machatts["compiler"] = compiler

        compsetnodes = self.get_nodes("compset")
        logger.debug("compsetnodes %d" % len(compsetnodes))
        for cnode in compsetnodes:
            gridnodes = self.get_nodes("grid", root=cnode)
            logger.debug("  gridnodes %d"%len(gridnodes))
            for gnode in gridnodes:
                testnamenodes = self.get_nodes("test",root=gnode)
                logger.debug("    testnamenodes %d"%len(testnamenodes))
                for tnode in testnamenodes:
                    machnodes = self.get_nodes("machine",machatts,root=tnode)
                    logger.debug("      machnodes %d"%len(machnodes))
                    for mach in machnodes:
                        thistest = {}
                        if machine is None or (machine is not None and mach.text == machine):
                            thistest["compiler"] = mach.get("compiler")
                            thistest["category"] = mach.get("testtype")
                            thistest["machine"] = mach.text
                            thistest["testname"] = tnode.get("name")
                            thistest["grid"] = gnode.get("name")
                            thistest["compset"] = cnode.get("name")
                            if ("testmods" in mach.attrib):
                                thistest["testmods"] = mach.get("testmods")
                            tests.append(thistest)
        return tests

    def _get_testsv2(self, machine=None, category=None, compiler=None):
        tests = []
        testnodes = self.get_nodes("test")
        machatts = {}
        if machine is not None:
            machatts["name"]      = machine
        if category is not None:
            machatts["category"]  = category
        if compiler is not None:
            machatts["compiler"]  = compiler

        for tnode in testnodes:
            machnodes = self.get_nodes("machine",machatts,root=tnode)
            if machnodes:
                this_test_node = {}
                for key, value in tnode.items():
                    if key == "name":
                        this_test_node["testname"] = value
                    else:
                        this_test_node[key] = value

                for mach in machnodes:
                    # this_test_node can include multiple tests
                    this_test = dict(this_test_node)
                    for key, value in mach.items():
                        if key == "name":
                            this_test["machine"] = value
                        else:
                            this_test[key] = value
                    # option xpath here gets only the children of tnode and ignores grandchildren
                    optionnodes = self.get_nodes("option", root=tnode, xpath="machines/options/option")

                    for onode in optionnodes:
                        this_test[onode.get('name')] = onode.text
                    optionnodes = self.get_nodes("option", root=mach)
                    for onode in optionnodes:
                        this_test[onode.get('name')] = onode.text


                    tests.append(this_test)

        return tests

    def get_tests(self, machine=None, category=None, compiler=None):
        if self.version == "1.0":
            return self._get_testsv1(machine, category, compiler)
        elif self.version == "2.0":
            return self._get_testsv2(machine, category, compiler)
        else:
            logger.critical("Did not recognize testlist file version %s for file %s"
                             % (self.version, self.filename))

