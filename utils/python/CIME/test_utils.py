# Utility functions used in create_test.py
from XML.standard_module_setup import *
from copy import deepcopy
from CIME.XML.testlist import Testlist

def  _get_tests_from_xml(xml_machine=None,xml_category=None,xml_compiler=None, xml_testlist=None,
                         machine=None, compiler=None):
    """
    Parse testlists for a list of tests
    """
    listoftests = []
    testlistfiles = []
    if(machine is not None):
        thismach=machine
    if(compiler is not None):
        thiscompiler = compiler

    if(xml_testlist is not None):
        expect(os.path.isfile(xml_testlist), "Testlist not found or not readable "+xml_testlist)
        testlistfiles.append(xml_testlist)
    else:
        files = Files()
        test_spec_files = files.get_values("TESTS_SPEC_FILE","component")
        for spec_file in test_spec_files.viewvalues():
            if(os.path.isfile(spec_file)):
                testlistfiles.append(spec_file)

    for testlistfile in testlistfiles:
        thistestlistfile = Testlist(testlistfile)
        newtests =  thistestlistfile.get_tests(xml_machine, xml_category, xml_compiler)
        for test in newtests:
            if(machine is None):
                thismach = test["machine"]
            if(compiler is None):
                thiscompiler = test["compiler"]
            test["name"] = "%s.%s.%s.%s_%s"%(test["testname"],test["grid"],test["compset"],thismach,thiscompiler)
            if ("testmods" in test):
                (moddir, modname) = test["testmods"].split("/")
                test["name"] += ".%s_%s"%(moddir, modname)
            logging.info("Adding test "+test["name"])
        listoftests += newtests

    return listoftests

def _convert_testlist_to_dict(test_names):
    from copy import deepcopy
    listoftests = []
    test = {}
    for name in test_names:
        test["name"] = name
        listoftests.append(deepcopy(test))
    return listoftests
