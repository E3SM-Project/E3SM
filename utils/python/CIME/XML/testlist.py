"""
Interface to the config_files.xml file.  This class inherits from generic_xml.py
It supports versions 1.0 and 2.0 of the testlist.xml file

In version 2 of the file options can be specified to further refine a test or
set of tests. They can be specified either at the top level, in which case they
apply to all machines/compilers for this test:

<test ...>
  <options>
    <option name="wallclock">00:20</option>
  </options>
  ...
</test>

or at the level of a particular machine/compiler:

<test ...>
  <machines>
    <machine ...>
      <options>
        <option name="wallclock">00:20</option>
      </options>
    </machine>
  </machines>
</test>

Currently supported options are:

- walltime: sets the wallclock limit in the queuing system

- memleak_tolerance: specifies the relative memory growth expected for this test

- comment: has no effect, but is written out when printing the test list

"""
from CIME.XML.standard_module_setup import *

from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files

logger = logging.getLogger(__name__)

class Testlist(GenericXML):

    def __init__(self,infile, files=None):
        """
        initialize an object
        """
        schema = None
        if files is None:
            files = Files()
        schema = files.get_schema("TESTS_SPEC_FILE")
        GenericXML.__init__(self, infile, schema=schema)

    def _get_testsv1(self, machine=None, category=None, compiler=None, compset=None, grid=None):
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
                        save_this = True
                        if machine is None or (machine is not None and mach.text == machine):
                            thistest["compiler"] = mach.get("compiler")
                            thistest["category"] = mach.get("testtype")
                            thistest["machine"] = mach.text
                            thistest["testname"] = tnode.get("name")
                            thistest["grid"] = gnode.get("name")
                            thistest["compset"] = cnode.get("name")
                            if ("testmods" in mach.attrib):
                                thistest["testmods"] = mach.get("testmods")
                            if compset is not None and compset != thistest["compset"]:
                                save_this = False
                            if grid is not None and grid != thistest["grid"]:
                                save_this = False
                            if save_this:
                                tests.append(thistest)
        return tests

    def _get_testsv2(self, machine=None, category=None, compiler=None, compset=None, grid=None):
        tests = []
        attributes = {}
        if compset is not None:
            attributes['compset'] = compset
        if grid is not None:
            attributes['grid'] = grid

        testnodes = self.get_nodes("test", attributes=attributes)

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
                    this_test["options"] = {}

                    # Get options that apply to all machines/compilers for this test
                    #
                    # option xpath here prevents us from looking within the
                    # machines block
                    optionnodes = self.get_nodes("option", root=tnode, xpath="options/option")
                    for onode in optionnodes:
                        this_test['options'][onode.get('name')] = onode.text

                    # Now get options specific to this machine/compiler
                    optionnodes = self.get_nodes("option", root=mach)
                    for onode in optionnodes:
                        this_test['options'][onode.get('name')] = onode.text


                    tests.append(this_test)

        return tests

    def get_tests(self, machine=None, category=None, compiler=None, compset=None, grid=None):
        if self.get_version() == 1.0:
            return self._get_testsv1(machine, category, compiler, compset, grid)
        elif self.get_version() >= 2.0:
            return self._get_testsv2(machine, category, compiler, compset, grid)
        else:
            logger.critical("Did not recognize testlist file version %s for file %s"
                             % (self.get_version(), self.filename))

