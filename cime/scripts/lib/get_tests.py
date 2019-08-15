import CIME.utils
from CIME.utils import expect, convert_to_seconds, parse_test_name, get_cime_root, get_model
from CIME.XML.machines import Machines
import six, sys, os

# Expect that, if a model wants to use python-based test lists, they will have a file
# config/$model/tests.py , containing a test dictionary called _TESTS

sys.path.insert(0, os.path.join(get_cime_root(), "config", get_model()))
_ALL_TESTS = {}
try:
    from tests import _TESTS # pylint: disable=import-error
    _ALL_TESTS.update(_TESTS)
except ImportError:
    pass

# Here are the tests belonging to cime suites. Format for individual tests is
# <test>.<grid>.<compset>[.<testmod>]
#
# suite_name : {
#     "inherit" : (suite1, suite2, ...), # Optional. Suites to inherit tests from. Default is None. Tuple, list, or str.
#     "time"    : "HH:MM:SS",            # Optional. Recommended upper-limit on test time.
#     "share"   : True|False,            # Optional. If True, all tests in this suite share a build. Default is False.
#     "tests"   : (test1, test2, ...)    # Optional. The list of tests for this suite. See above for format. Tuple, list, or str. This is the ONLY inheritable attribute.
# }

_CIME_TESTS = {

    "cime_tiny" : {
        "time"  : "0:10:00",
        "tests" : (
            "ERS.f19_g16_rx1.A",
            "NCK.f19_g16_rx1.A",
            )
        },

    "cime_test_only_pass" : {
        "time"  : "0:10:00",
        "tests" : (
            "TESTRUNPASS_P1.f19_g16_rx1.A",
            "TESTRUNPASS_P1.ne30_g16_rx1.A",
            "TESTRUNPASS_P1.f45_g37_rx1.A",
            )
        },

    "cime_test_only_slow_pass" : {
        "time"  : "0:10:00",
        "tests" : (
            "TESTRUNSLOWPASS_P1.f19_g16_rx1.A",
            "TESTRUNSLOWPASS_P1.ne30_g16_rx1.A",
            "TESTRUNSLOWPASS_P1.f45_g37_rx1.A",
            )
        },

    "cime_test_only" : {
        "time"  : "0:10:00",
        "tests" : (
            "TESTBUILDFAIL_P1.f19_g16_rx1.A",
            "TESTBUILDFAILEXC_P1.f19_g16_rx1.A",
            "TESTRUNFAIL_P1.f19_g16_rx1.A",
            "TESTRUNSTARCFAIL_P1.f19_g16_rx1.A",
            "TESTRUNFAILEXC_P1.f19_g16_rx1.A",
            "TESTRUNPASS_P1.f19_g16_rx1.A",
            "TESTTESTDIFF_P1.f19_g16_rx1.A",
            "TESTMEMLEAKFAIL_P1.f09_g16.X",
            "TESTMEMLEAKPASS_P1.f09_g16.X",
            )
        },

    "cime_test_all" : {
        "inherit" : "cime_test_only",
        "time"    : "0:10:00",
        "tests"   : "TESTRUNDIFF_P1.f19_g16_rx1.A"
        },

    "cime_test_share" : {
        "time"  : "0:10:00",
        "share" : True,
        "tests" : (
            "SMS_P2.f19_g16_rx1.A",
            "SMS_P4.f19_g16_rx1.A",
            "SMS_P8.f19_g16_rx1.A",
            "SMS_P16.f19_g16_rx1.A",
            )
        },

    "cime_test_share2" : {
        "time"  : "0:10:00",
        "share" : True,
        "tests" : (
            "SMS_P2.f19_g16_rx1.X",
            "SMS_P4.f19_g16_rx1.X",
            "SMS_P8.f19_g16_rx1.X",
            "SMS_P16.f19_g16_rx1.X",
            )
        },

    "cime_test_repeat" : {
        "tests" : (
            "TESTRUNPASS_P1.f19_g16_rx1.A",
            "TESTRUNPASS_P2.ne30_g16_rx1.A",
            "TESTRUNPASS_P4.f45_g37_rx1.A",
            )
        },

    "cime_test_time" : {
        "time"    : "0:13:00",
        "tests" : (
            "TESTRUNPASS_P69.f19_g16_rx1.A.testmod",
            )
        },

    "cime_test_multi_inherit" : {
        "inherit" : ("cime_test_repeat", "cime_test_only_pass", "cime_test_all")
        },

    "cime_developer" : {
        "time"  : "0:15:00",
        "tests" : (
            "NCK_Ld3.f45_g37_rx1.A",
            "ERI.f09_g16.X",
            "ERIO.f09_g16.X",
            "SEQ_Ln9.f19_g16_rx1.A",
            "ERS.ne30_g16_rx1.A.drv-y100k",
            "IRT_N2.f19_g16_rx1.A",
            "ERR.f45_g37_rx1.A",
            "ERP.f45_g37_rx1.A",
            "SMS_D_Ln9_Mmpi-serial.f19_g16_rx1.A",
            "DAE.ww3a.ADWAV",
            "PET_P4.f19_f19.A",
            "PEM_P4.f19_f19.A",
            "SMS.T42_T42.S",
            "PRE.f19_f19.ADESP",
            "PRE.f19_f19.ADESP_TEST",
            "MCC_P1.f19_g16_rx1.A",
            "LDSTA.f45_g37_rx1.A",
            )
        },

}

_ALL_TESTS.update(_CIME_TESTS)

###############################################################################
def _get_key_data(raw_dict, key, the_type):
###############################################################################
    if key not in raw_dict:
        if the_type is tuple:
            return ()
        elif the_type is str:
            return None
        elif the_type is bool:
            return False
        else:
            expect(False, "Unsupported type {}".format(the_type))
    else:
        val = raw_dict[key]
        if the_type is tuple and isinstance(val, six.string_types):
            val = (val, )

        expect(isinstance(val, the_type),
               "Wrong type for {}, {} is a {} but expected {}".format(key, val, type(val), the_type))

        return val

###############################################################################
def get_test_data(suite):
###############################################################################
    """
    For a given suite, returns (inherit, time, share, tests)
    """
    raw_dict = _ALL_TESTS[suite]
    for key in raw_dict.keys():
        expect(key in ["inherit", "time", "share", "tests"], "Unexpected test key '{}'".format(key))

    return _get_key_data(raw_dict, "inherit", tuple), _get_key_data(raw_dict, "time", str), _get_key_data(raw_dict, "share", bool), _get_key_data(raw_dict, "tests", tuple)

###############################################################################
def get_test_suites():
###############################################################################
    return list(_ALL_TESTS.keys())

###############################################################################
def get_test_suite(suite, machine=None, compiler=None, skip_inherit=False):
###############################################################################
    """
    Return a list of FULL test names for a suite.
    """
    expect(suite in get_test_suites(), "Unknown test suite: '{}'".format(suite))
    machobj = Machines(machine=machine)
    machine = machobj.get_machine_name()

    if(compiler is None):
        compiler = machobj.get_default_compiler()
    expect(machobj.is_valid_compiler(compiler),"Compiler {} not valid for machine {}".format(compiler,machine))

    inherits_from, _, _, tests_raw = get_test_data(suite)
    tests = []
    for item in tests_raw:
        expect(isinstance(item, six.string_types), "Bad type of test {}, expected string".format(item))

        test_mod = None
        test_components = item.split(".")
        expect(len(test_components) in [3, 4], "Bad test name {}".format(item))

        if (len(test_components) == 4):
            test_name = ".".join(test_components[:-1])
            test_mod = test_components[-1]
        else:
            test_name = item

        tests.append(CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler, testmod=test_mod))

    if not skip_inherit:
        for inherits in inherits_from:
            inherited_tests = get_test_suite(inherits, machine, compiler)

            for inherited_test in inherited_tests:
                if inherited_test not in tests:
                    tests.append(inherited_test)

    return tests

###############################################################################
def suite_has_test(suite, test_full_name, skip_inherit=False):
###############################################################################
    _, _, _, _, machine, compiler, _ = CIME.utils.parse_test_name(test_full_name)
    expect(machine is not None, "{} is not a full test name".format(test_full_name))

    tests = get_test_suite(suite, machine=machine, compiler=compiler, skip_inherit=skip_inherit)
    return test_full_name in tests

###############################################################################
def get_build_groups(tests):
###############################################################################
    """
    Given a list of tests, return a list of lists, with each list representing
    a group of tests that can share executables.

    >>> tests = ["SMS_P2.f19_g16_rx1.A.melvin_gnu", "SMS_P4.f19_g16_rx1.A.melvin_gnu", "SMS_P2.f19_g16_rx1.X.melvin_gnu", "SMS_P4.f19_g16_rx1.X.melvin_gnu", "TESTRUNSLOWPASS_P1.f19_g16_rx1.A.melvin_gnu", "TESTRUNSLOWPASS_P1.ne30_g16_rx1.A.melvin_gnu"]
    >>> get_build_groups(tests)
    [('SMS_P2.f19_g16_rx1.A.melvin_gnu', 'SMS_P4.f19_g16_rx1.A.melvin_gnu'), ('SMS_P2.f19_g16_rx1.X.melvin_gnu', 'SMS_P4.f19_g16_rx1.X.melvin_gnu'), ('TESTRUNSLOWPASS_P1.f19_g16_rx1.A.melvin_gnu',), ('TESTRUNSLOWPASS_P1.ne30_g16_rx1.A.melvin_gnu',)]
    """
    build_groups = [] # list of tuples ([tests], set(suites))

    # Get a list of suites that share exes
    suites = get_test_suites()
    share_suites = []
    for suite in suites:
        share = get_test_data(suite)[2]
        if share:
            share_suites.append(suite)

    # Divide tests up into build groups. Assumes that build-compatibility is transitive
    for test in tests:
        matched = False

        my_share_suites = set()
        for suite in share_suites:
            if suite_has_test(suite, test, skip_inherit=True):
                my_share_suites.add(suite)

        # Try to match this test with an existing build group
        if my_share_suites:
            for build_group_tests, build_group_suites in build_groups:
                overlap = build_group_suites & my_share_suites
                if overlap:
                    matched = True
                    build_group_tests.append(test)
                    build_group_suites.update(my_share_suites)
                    break

        # Nothing matched, this test is in a build group of its own
        if not matched:
            build_groups.append(([test], my_share_suites))

    return [tuple(item[0]) for item in build_groups]

###############################################################################
def infer_machine_name_from_tests(testargs):
###############################################################################
    """
    >>> infer_machine_name_from_tests(["NCK.f19_g16_rx1.A.melvin_gnu"])
    'melvin'
    >>> infer_machine_name_from_tests(["NCK.f19_g16_rx1.A"])
    >>> infer_machine_name_from_tests(["NCK.f19_g16_rx1.A", "NCK.f19_g16_rx1.A.melvin_gnu"])
    'melvin'
    >>> infer_machine_name_from_tests(["NCK.f19_g16_rx1.A.melvin_gnu", "NCK.f19_g16_rx1.A.melvin_gnu"])
    'melvin'
    """
    e3sm_test_suites = get_test_suites()

    machine = None
    for testarg in testargs:
        testarg = testarg.strip()
        if testarg.startswith("^"):
            testarg = testarg[1:]

        if testarg not in e3sm_test_suites:
            machine_for_this_test = parse_test_name(testarg)[4]
            if machine_for_this_test is not None:
                if machine is None:
                    machine = machine_for_this_test
                else:
                    expect(machine == machine_for_this_test, "Must have consistent machine '%s' != '%s'" % (machine, machine_for_this_test))

    return machine

###############################################################################
def get_full_test_names(testargs, machine, compiler):
###############################################################################
    """
    Return full test names in the form:
    TESTCASE.GRID.COMPSET.MACHINE_COMPILER.TESTMODS
    Testmods are optional

    Testargs can be categories or test names and support the NOT symbol '^'

    >>> get_full_test_names(["cime_tiny"], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu', 'NCK.f19_g16_rx1.A.melvin_gnu']

    >>> get_full_test_names(["cime_tiny", "PEA_P1_M.f45_g37_rx1.A"], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu', 'NCK.f19_g16_rx1.A.melvin_gnu', 'PEA_P1_M.f45_g37_rx1.A.melvin_gnu']

    >>> get_full_test_names(['ERS.f19_g16_rx1.A', 'NCK.f19_g16_rx1.A', 'PEA_P1_M.f45_g37_rx1.A'], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu', 'NCK.f19_g16_rx1.A.melvin_gnu', 'PEA_P1_M.f45_g37_rx1.A.melvin_gnu']

    >>> get_full_test_names(["cime_tiny", "^NCK.f19_g16_rx1.A"], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu']

    >>> get_full_test_names(["cime_test_multi_inherit"], "melvin", "gnu")
    ['TESTBUILDFAILEXC_P1.f19_g16_rx1.A.melvin_gnu', 'TESTBUILDFAIL_P1.f19_g16_rx1.A.melvin_gnu', 'TESTMEMLEAKFAIL_P1.f09_g16.X.melvin_gnu', 'TESTMEMLEAKPASS_P1.f09_g16.X.melvin_gnu', 'TESTRUNDIFF_P1.f19_g16_rx1.A.melvin_gnu', 'TESTRUNFAILEXC_P1.f19_g16_rx1.A.melvin_gnu', 'TESTRUNFAIL_P1.f19_g16_rx1.A.melvin_gnu', 'TESTRUNPASS_P1.f19_g16_rx1.A.melvin_gnu', 'TESTRUNPASS_P1.f45_g37_rx1.A.melvin_gnu', 'TESTRUNPASS_P1.ne30_g16_rx1.A.melvin_gnu', 'TESTRUNPASS_P2.ne30_g16_rx1.A.melvin_gnu', 'TESTRUNPASS_P4.f45_g37_rx1.A.melvin_gnu', 'TESTRUNSTARCFAIL_P1.f19_g16_rx1.A.melvin_gnu', 'TESTTESTDIFF_P1.f19_g16_rx1.A.melvin_gnu']
    """
    expect(machine is not None, "Must define a machine")
    expect(compiler is not None, "Must define a compiler")
    e3sm_test_suites = get_test_suites()

    tests_to_run = set()
    negations = set()

    for testarg in testargs:
        # remove any whitespace in name
        testarg = testarg.strip()
        if (testarg.startswith("^")):
            negations.add(testarg[1:])
        elif (testarg in e3sm_test_suites):
            tests_to_run.update(get_test_suite(testarg, machine, compiler))
        else:
            try:
                tests_to_run.add(CIME.utils.get_full_test_name(testarg, machine=machine, compiler=compiler))
            except Exception:
                if "." not in testarg:
                    expect(False, "Unrecognized test suite '{}'".format(testarg))
                else:
                    raise

    for negation in negations:
        if (negation in e3sm_test_suites):
            tests_to_run -= set(get_test_suite(negation, machine, compiler))
        else:
            fullname = CIME.utils.get_full_test_name(negation, machine=machine, compiler=compiler)
            if (fullname in tests_to_run):
                tests_to_run.remove(fullname)

    return list(sorted(tests_to_run))

###############################################################################
def get_recommended_test_time(test_full_name):
###############################################################################
    """
    >>> get_recommended_test_time("ERS.f19_g16_rx1.A.melvin_gnu")
    '0:10:00'

    >>> get_recommended_test_time("TESTRUNPASS_P69.f19_g16_rx1.A.melvin_gnu.testmod")
    '0:13:00'

    >>> get_recommended_test_time("PET_Ln20.ne30_ne30.FC5.sandiatoss3_intel.cam-outfrq9s")
    >>>
    """
    best_time = None
    suites = get_test_suites()
    for suite in suites:
        rec_time = get_test_data(suite)[1]
        if suite_has_test(suite, test_full_name, skip_inherit=True) and rec_time is not None and \
           (best_time is None or convert_to_seconds(rec_time) < convert_to_seconds(best_time)):
            best_time = rec_time

    return best_time

###############################################################################
def key_test_time(test_full_name):
###############################################################################
    result = get_recommended_test_time(test_full_name)
    return 99999999 if result is None else convert_to_seconds(result)
