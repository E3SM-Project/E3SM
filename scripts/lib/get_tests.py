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
except:
    pass

# Here are the tests belonging to e3sm suites. Format is
# <test>.<grid>.<compset>.
# suite_name -> (inherits_from, timelimit, [test [, mods[, machines]]])
#   To elaborate, if no mods are needed, a string representing the testname is all that is needed.
#   If testmods are needed, a 2-ple must be provided  (test, mods)
#   If you want to restrict the test mods to certain machines, than a 3-ple is needed (test, mods, [machines])

_CIME_TESTS = {

    "cime_tiny" : (None, "0:10:00",
                   ("ERS.f19_g16_rx1.A",
                    "NCK.f19_g16_rx1.A")
                   ),

    "cime_test_only_pass" : (None, "0:10:00",
                   ("TESTRUNPASS_P1.f19_g16_rx1.A",
                    "TESTRUNPASS_P1.ne30_g16_rx1.A",
                    "TESTRUNPASS_P1.f45_g37_rx1.A")
                   ),

    "cime_test_only_slow_pass" : (None, "0:10:00",
                   ("TESTRUNSLOWPASS_P1.f19_g16_rx1.A",
                    "TESTRUNSLOWPASS_P1.ne30_g16_rx1.A",
                    "TESTRUNSLOWPASS_P1.f45_g37_rx1.A")
                   ),

    "cime_test_only" : (None, "0:10:00",
                   ("TESTBUILDFAIL_P1.f19_g16_rx1.A",
                    "TESTBUILDFAILEXC_P1.f19_g16_rx1.A",
                    "TESTRUNFAIL_P1.f19_g16_rx1.A",
                    "TESTRUNSTARCFAIL_P1.f19_g16_rx1.A",
                    "TESTRUNFAILEXC_P1.f19_g16_rx1.A",
                    "TESTRUNPASS_P1.f19_g16_rx1.A",
                    "TESTTESTDIFF_P1.f19_g16_rx1.A",
                    "TESTMEMLEAKFAIL_P1.f09_g16.X",
                    "TESTMEMLEAKPASS_P1.f09_g16.X")
                   ),

    "cime_test_all" : ("cime_test_only", "0:10:00",
                       ("TESTRUNDIFF_P1.f19_g16_rx1.A", )),

    "cime_developer" : (None, "0:15:00",
                            ("NCK_Ld3.f45_g37_rx1.A",
                             "ERI.f09_g16.X",
                             "ERIO.f09_g16.X",
                             "SEQ_Ln9.f19_g16_rx1.A",
                             ("ERS.ne30_g16_rx1.A","drv-y100k"),
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
                             "LDSTA.f45_g37_rx1.A")
                            ),

}

_ALL_TESTS.update(_CIME_TESTS)

###############################################################################
def get_test_suite(suite, machine=None, compiler=None):
###############################################################################
    """
    Return a list of FULL test names for a suite.
    """
    expect(suite in _ALL_TESTS, "Unknown test suite: '{}'".format(suite))
    machobj = Machines(machine=machine)
    machine = machobj.get_machine_name()

    if(compiler is None):
        compiler = machobj.get_default_compiler()
    expect(machobj.is_valid_compiler(compiler),"Compiler {} not valid for machine {}".format(compiler,machine))

    inherits_from, _, tests_raw = _ALL_TESTS[suite]
    tests = []
    for item in tests_raw:
        test_mod = None
        if (isinstance(item, six.string_types)):
            test_name = item
        else:
            expect(isinstance(item, tuple), "Bad item type for item '{}'".format(str(item)))
            expect(len(item) in [2, 3], "Expected two or three items in item '{}'".format(str(item)))
            expect(isinstance(item[0], six.string_types), "Expected string in first field of item '{}'".format(str(item)))
            expect(isinstance(item[1], six.string_types), "Expected string in second field of item '{}'".format(str(item)))

            test_name = item[0]
            if (len(item) == 2):
                test_mod = item[1]
            else:
                expect(type(item[2]) in [six.string_types, tuple], "Expected string or tuple for third field of item '{}'".format(str(item)))
                test_mod_machines = [item[2]] if isinstance(item[2], six.string_types) else item[2]
                if (machine in test_mod_machines):
                    test_mod = item[1]

        tests.append(CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler, testmod=test_mod))

    if (inherits_from is not None):
        inherits_from = [inherits_from] if isinstance(inherits_from, six.string_types) else inherits_from
        for inherits in inherits_from:
            inherited_tests = get_test_suite(inherits, machine, compiler)

            expect(len(set(tests) & set(inherited_tests)) == 0,
                   "Tests {} defined in multiple suites".format(", ".join(set(tests) & set(inherited_tests))))
            tests.extend(inherited_tests)

    return tests

###############################################################################
def get_test_suites():
###############################################################################
    return list(_ALL_TESTS.keys())

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
            except:
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

    >>> get_recommended_test_time("ERP_Ln9.ne30_ne30.FC5.melvin_gun.cam-outfrq9s")
    '0:45:00'

    >>> get_recommended_test_time("PET_Ln9.ne30_ne30.FC5.sandiatoss3_intel.cam-outfrq9s")
    '03:00:00'

    >>> get_recommended_test_time("PET_Ln20.ne30_ne30.FC5.sandiatoss3_intel.cam-outfrq9s")
    >>>
    """
    _, _, _, _, machine, compiler, _ = CIME.utils.parse_test_name(test_full_name)
    expect(machine is not None, "{} is not a full test name".format(test_full_name))

    best_time = None
    suites = get_test_suites()
    for suite in suites:
        _, rec_time, tests_raw = _ALL_TESTS[suite]
        for item in tests_raw:
            test_mod = None
            if (isinstance(item, six.string_types)):
                test_name = item
            else:
                test_name = item[0]
                if (len(item) == 2):
                    test_mod = item[1]
                else:
                    test_mod_machines = [item[2]] if isinstance(item[2], six.string_types) else item[2]
                    if (machine in test_mod_machines):
                        test_mod = item[1]

            full_test = CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler, testmod=test_mod)

            if full_test == test_full_name and rec_time is not None:
                if best_time is None or \
                        convert_to_seconds(rec_time) < convert_to_seconds(best_time):
                    best_time = rec_time

    return best_time

###############################################################################
def sort_by_time(test_one, test_two):
###############################################################################
    """
    >>> tests = get_full_test_names(["cime_tiny"], "melvin", "gnu")
    >>> tests.extend(get_full_test_names(["cime_developer"], "melvin", "gnu"))
    >>> tests.append("A.f19_f19.A.melvin_gnu")
    >>> tests.sort(cmp=sort_by_time)
    >>> tests
    ['DAE.f19_f19.A.melvin_gnu', 'ERI.f09_g16.X.melvin_gnu', 'ERIO.f09_g16.X.melvin_gnu', 'ERP.f45_g37_rx1.A.melvin_gnu', 'ERR.f45_g37_rx1.A.melvin_gnu', 'ERS.ne30_g16_rx1.A.melvin_gnu', 'IRT_N2.f19_g16_rx1.A.melvin_gnu', 'NCK_Ld3.f45_g37_rx1.A.melvin_gnu', 'PET_P32.f19_f19.A.melvin_gnu', 'PRE.f19_f19.ADESP.melvin_gnu', 'PRE.f19_f19.ADESP_TEST.melvin_gnu', 'SEQ_Ln9.f19_g16_rx1.A.melvin_gnu', 'SMS.T42_T42.S.melvin_gnu', 'SMS_D_Ln9.f19_g16_rx1.A.melvin_gnu', 'ERS.f19_g16_rx1.A.melvin_gnu', 'NCK.f19_g16_rx1.A.melvin_gnu', 'A.f19_f19.A.melvin_gnu']
    """
    rec1, rec2 = get_recommended_test_time(test_one), get_recommended_test_time(test_two)
    if rec1 == rec2:
        return (test_one > test_two) - (test_two < test_one)
    else:
        if rec2 is None:
            return -1
        elif rec1 is None:
            return 1
        else:
            a = convert_to_seconds(rec2)
            b = convert_to_seconds(rec1)
            return (a < b) - (b < a)
