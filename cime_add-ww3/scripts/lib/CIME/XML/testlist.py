"""
Interface to the config_files.xml file.  This class inherits from generic_xml.py
It supports version 2.0 of the testlist.xml file

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

- workflow: adds a workflow to the test
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
        expect(self.get_version() >= 2.0,
               "{} is an unsupported version of the testfile format and will be ignored".format(infile))

    def get_tests(self, machine=None, category=None, compiler=None, compset=None, grid=None, supported_only=False):
        tests = []
        attributes = {}
        if compset is not None:
            attributes['compset'] = compset
        if grid is not None:
            attributes['grid'] = grid

        testnodes = self.get_children("test", attributes=attributes)

        machatts = {}
        if machine is not None:
            machatts["name"]      = machine
        if category is not None:
            machatts["category"]  = category
        if compiler is not None:
            machatts["compiler"]  = compiler


        for tnode in testnodes:
            if supported_only and self.has(tnode, "supported") and self.get(tnode, "supported") == 'false':
                continue

            machnode = self.get_optional_child("machines", root=tnode)
            machnodes = None if machnode is None else self.get_children("machine",machatts,root=machnode)
            if machnodes:
                this_test_node = {}
                for key, value in self.attrib(tnode).items():
                    if key == "name":
                        this_test_node["testname"] = value
                    else:
                        this_test_node[key] = value



                # Get options that apply to all machines/compilers for this test
                options = self.get_children("options", root=tnode)
                if len(options) > 0:
                    optionnodes = self.get_children("option", root=options[0])
                else:
                    optionnodes = []
                for mach in machnodes:
                    # this_test_node can include multiple tests
                    this_test = dict(this_test_node)
                    for key, value in self.attrib(mach).items():
                        if key == "name":
                            this_test["machine"] = value
                        else:
                            this_test[key] = value
                    this_test["options"] = {}

                    for onode in optionnodes:
                        this_test['options'][self.get(onode, 'name')] = self.text(onode)

                    # Now get options specific to this machine/compiler
                    options = self.get_optional_child("options", root=mach)
                    optionnodes = [] if options is None else self.get_children("option", root=options)
                    for onode in optionnodes:
                        this_test['options'][self.get(onode, 'name')] = self.text(onode)

                    tests.append(this_test)

        return tests
