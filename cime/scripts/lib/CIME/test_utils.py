"""
Utility functions used in test_scheduler.py, and by other utilities that need to
get test lists.
"""
import glob
from CIME.XML.standard_module_setup import *
from CIME.XML.testlist import Testlist
from CIME.XML.files import Files
from CIME.test_status import TEST_STATUS_FILENAME
import CIME.utils

logger = logging.getLogger(__name__)

def get_tests_from_xml(xml_machine=None,xml_category=None,xml_compiler=None, xml_testlist=None,
                       machine=None, compiler=None, driver=None):
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
        comps = files.get_components("TESTS_SPEC_FILE")
        for comp in comps:
            test_spec_file = files.get_value("TESTS_SPEC_FILE", {"component":comp})
            if(os.path.isfile(test_spec_file)):
                testlistfiles.append(test_spec_file)

    for testlistfile in testlistfiles:
        thistestlistfile = Testlist(testlistfile)
        logger.debug("Testlist file is "+testlistfile)
        logger.debug("xml_machine {} xml_category {} xml_compiler {}".format(xml_machine, xml_category, xml_compiler))
        newtests =  thistestlistfile.get_tests(xml_machine, xml_category, xml_compiler)
        for test in newtests:
            if(machine is None):
                thismach = test["machine"]
            if(compiler is None):
                thiscompiler = test["compiler"]
            test["name"] = CIME.utils.get_full_test_name(test["testname"], grid=test["grid"], compset=test["compset"],
                                                         machine=thismach, compiler=thiscompiler,
                                                         testmod=None if "testmods" not in test else test["testmods"])
            if driver:
                # override default or specified driver
                founddriver = False
                for specdriver in ("Vnuopc","Vmct","Vmoab"):
                    if specdriver in test["name"]:
                        test["name"] = test["name"].replace(specdriver,"V{}".format(driver))
                        founddriver = True
                if not founddriver:
                    name = test["name"]
                    index = name.find('.')
                    test["name"] = name[:index] + "_V{}".format(driver) + name[index:]


            logger.debug("Adding test {} with compiler {}".format(test["name"], test["compiler"]))
        listoftests += newtests
        logger.debug("Found {:d} tests".format(len(listoftests)))

    return listoftests

def test_to_string(test, category_field_width=0, test_field_width=0, show_options=False):
    """Given a test dictionary, return a string representation suitable for printing

    Args:
        test (dict): dictionary for a single test - e.g., one element from the
                     list returned by get_tests_from_xml
        category_field_width (int): minimum amount of space to use for printing the test category
        test_field_width (int): minimum amount of space to use for printing the test category
        show_options (bool): if True, print test options, too (note that the 'comment'
                             option is always printed, if present)

    Basic functionality:
    >>> mytest = {'name': 'SMS.f19_g16.A.cheyenne_intel', 'category': 'prealpha', 'options': {}}
    >>> test_to_string(mytest, 10)
    'prealpha  : SMS.f19_g16.A.cheyenne_intel'

    Printing comments:
    >>> mytest = {'name': 'SMS.f19_g16.A.cheyenne_intel', 'category': 'prealpha', 'options': {'comment': 'my remarks'}}
    >>> test_to_string(mytest, 10)
    'prealpha  : SMS.f19_g16.A.cheyenne_intel # my remarks'

    Printing other options, too:
    >>> mytest = {'name': 'SMS.f19_g16.A.cheyenne_intel', 'category': 'prealpha', 'options': {'comment': 'my remarks', 'wallclock': '0:20', 'memleak_tolerance': 0.2}}
    >>> test_to_string(mytest, 10, show_options=True)
    'prealpha  : SMS.f19_g16.A.cheyenne_intel # my remarks  # memleak_tolerance: 0.2  # wallclock: 0:20'
    """

    mystr = "%-*s: %-*s"%(category_field_width, test['category'], test_field_width, test['name'])
    if 'options' in test:
        myopts = test['options'].copy()
        comment = myopts.pop('comment', None)
        if comment:
            mystr += " # {}".format(comment)
        if show_options:
            for one_opt in sorted(myopts):
                mystr += "  # {}: {}".format(one_opt, myopts[one_opt])

    return mystr

def get_test_status_files(test_root, compiler, test_id=None):
    test_id_glob = "*{}*".format(compiler) if test_id is None else "*{}*".format(test_id)
    test_status_files = glob.glob("{}/{}/{}".format(test_root, test_id_glob, TEST_STATUS_FILENAME))
    test_status_files = [item for item in test_status_files if not os.path.dirname(item).endswith('ref1') and not os.path.dirname(item).endswith('ref2')]

    expect(test_status_files, "No matching test cases found in for {}/{}/{}".format(test_root, test_id_glob, TEST_STATUS_FILENAME))
    return test_status_files
