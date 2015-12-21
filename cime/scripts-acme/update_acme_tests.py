import os, tempfile

import acme_util
from acme_util import expect, warning

# Here are the tests belonging to acme suites. Format is
# <test>.<grid>.<compset>.
# suite_name -> (inherits_from, [test [, mods[, machines]]])
#   To elaborate, if no mods are needed, a string representing the testname is all that is needed.
#   If testmods are needed, a 2-ple must be provided  (test, mods)
#   If you want to restrict the test mods to certain machines, than a 3-ple is needed (test, mods, [machines])
_TEST_SUITES = {
    "acme_tiny" : (None,
                   ("ERS.f19_g16_rx1.A",
                    "NCK.f19_g16_rx1.A")
                   ),

    "acme_test_only_pass" : (None,
                   ("TESTRUNPASS_P1.f19_g16_rx1.A",
                    "TESTRUNPASS_P1.ne30_g16_rx1.A",
                    "TESTRUNPASS_P1.f45_g37_rx1.A")
                   ),

    "acme_test_only_slow_pass" : (None,
                   ("TESTRUNSLOWPASS_P1.f19_g16_rx1.A",
                    "TESTRUNSLOWPASS_P1.ne30_g16_rx1.A",
                    "TESTRUNSLOWPASS_P1.f45_g37_rx1.A")
                   ),

    "acme_test_only" : (None,
                   ("TESTBUILDFAIL.f19_g16_rx1.A",
                    "TESTRUNFAIL_P1.f19_g16_rx1.A",
                    "TESTRUNPASS_P1.f19_g16_rx1.A")
                   ),

    "acme_land_developer" : (None,
                             ("SMS.f19_f19.I1850CLM45CN",
                              "SMS.f09_g16.I1850CLM45CN",
                              "SMS.hcru_hcru.I1850CRUCLM45CN")
                             ),

    "acme_atm_developer" : (None,
                            ("ERS.ne16_ne16.FC5MAM4",
                             "ERS.ne16_ne16.FC5PLMOD",
                             "ERS.ne16_ne16.FC5CLBMG2",
                             "ERS.ne16_ne16.FC5CLBMG2MAM4",
                             "ERS.ne16_ne16.FC5CLBMG2MAM4RESUS",
                             "ERS.ne16_ne16.FC5CLBMG2MAM4RESUSCOSP",
                             "ERS.ne16_ne16.FC5CLBMG2MAM4RESUSCOSPBC",
                             "ERS.ne16_ne16.FC5ATMMOD")
                            ),

    "acme_developer" : ("acme_land_developer",
                        ("ERS.f19_g16_rx1.A",
                         "ERS.f45_g37_rx1.DTEST",
                         "ERS.ne30_g16_rx1.A",
                         "ERS_IOP.f19_g16_rx1.A",
                         "ERS_IOP.f45_g37_rx1.DTEST",
                         "ERS_IOP.ne30_g16_rx1.A",
                         "ERS_IOP4c.f19_g16_rx1.A",
                         "ERS_IOP4c.ne30_g16_rx1.A",
                         "ERS_IOP4p.f19_g16_rx1.A",
                         "ERS_IOP4p.ne30_g16_rx1.A",
                         ("ERP_Ln9.ne30_ne30.FC5", "cam-outfrq9s"),
                         "NCK.f19_g16_rx1.A",
                         "PEA_P1_M.f45_g37_rx1.A",
                         "SMS.ne30_f19_g16_rx1.A",
                         "ERS_Ld5.T62_mpas120.C_MPAS_NORMAL_YEAR",
                         "ERS.f09_g16_g.MPASLI_ONLY",
                         "SMS.T62_mpas120_gis20.MPAS_LISIO_TEST",
                         "SMS.f09_g16_a.IGCLM45_MLI")
                        ),

    "acme_integration" : ("acme_developer",
                          ("ERS.ne30_g16.B1850C5",
                           "ERS.f19_f19.FAMIPC5",
                           "ERS.ne16_ne16.FC5PM",
                           "ERS.ne16_g37.B1850C5",
                           "ERS.ne16_ne16.FC5PLMOD",
                           "ERS.ne16_ne16.FC5MAM4",
                           "ERS.f45_g37.B1850C5",
                           "ERS_D.f45_g37.B1850C5",
                           "ERS_IOP_Ld3.f19_f19.FAMIPC5",
                           "ERS_Ld3.ne16_g37.FC5",
                           "ERS_Ld3.ne30_ne30.FC5",
                           "ERT.ne16_g37.B1850C5",
                           ("PET_PT_Ln9.ne30_ne30.FC5", "cam-outfrq9s"),
                           "PET_PT.f19_g16.X",
                           "PET_PT.f45_g37_rx1.A",
                           "PFS.ne30_ne30.FC5",
                           "SEQ_IOP_PFC.f19_g16.X",
                           "SEQ_PFC.f45_g37.B1850C5",
                           "SMS.ne30_m120.A_B1850CN",
                           "SMS.ne16_ne16.FC5AQUAP",
                           "SMS_D.f19_g16.B20TRC5",
                           "SMS_D_Ld3.ne16_ne16.FC5",
                           "SMS.f09_g16_a.MPASLIALB_ONLY")
                          ),
}

###############################################################################
def get_test_suite(suite, machine=None, compiler=None):
###############################################################################
    """
    Return a list of FULL test names for a suite.
    """
    expect(suite in _TEST_SUITES, "Unknown test suite: '%s'" % suite)

    machine = acme_util.probe_machine_name() if machine is None else machine
    compiler = acme_util.get_machine_info("COMPILERS", machine=machine)[0] if compiler is None else compiler

    inherits_from, tests_raw = _TEST_SUITES[suite]
    tests = []
    for item in tests_raw:
        test_mod = None
        if (isinstance(item, str)):
            test_name = item
        else:
            expect(isinstance(item, tuple), "Bad item type for item '%s'" % str(item))
            expect(len(item) in [2, 3], "Expected two or three items in item '%s'" % str(item))
            expect(isinstance(item[0], str), "Expected string in first field of item '%s'" % str(item))
            expect(isinstance(item[1], str), "Expected string in second field of item '%s'" % str(item))

            test_name = item[0]
            if (len(item) == 2):
                test_mod = item[1]
            else:
                expect(type(item[2]) in [str, tuple], "Expected string or tuple for third field of item '%s'" % str(item))
                test_mod_machines = [item[2]] if isinstance(item[2], str) else item[2]
                if (machine in test_mod_machines):
                    test_mod = item[1]

        tests.append(acme_util.get_full_test_name(test_name, machine, compiler, testmod=test_mod))

    if (inherits_from is not None):
        inherited_tests = get_test_suite(inherits_from, machine, compiler)

        expect(len(set(tests) & set(inherited_tests)) == 0,
               "Tests %s defined in multiple suites" % ", ".join(set(tests) & set(inherited_tests)))
        tests.extend(inherited_tests)

    return tests

###############################################################################
def get_test_suites():
###############################################################################
    return _TEST_SUITES.keys()

###############################################################################
def get_full_test_names(testargs, machine, compiler):
###############################################################################
    """
    Return full test names in the form:
    TESTCASE.GRID.COMPSET.MACHINE_COMPILER.TESTMODS
    Testmods are optional

    Testargs can be categories or test names and support the NOT symbol '^'

    >>> get_full_test_names(["acme_tiny"], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu', 'NCK.f19_g16_rx1.A.melvin_gnu']

    >>> get_full_test_names(["acme_tiny"], "melvin", "intel")
    ['ERS.f19_g16_rx1.A.melvin_intel', 'NCK.f19_g16_rx1.A.melvin_intel']

    >>> get_full_test_names(["acme_tiny", "PEA_P1_M.f45_g37_rx1.A"], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu', 'NCK.f19_g16_rx1.A.melvin_gnu', 'PEA_P1_M.f45_g37_rx1.A.melvin_gnu']

    >>> get_full_test_names(['ERS.f19_g16_rx1.A', 'NCK.f19_g16_rx1.A', 'PEA_P1_M.f45_g37_rx1.A'], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu', 'NCK.f19_g16_rx1.A.melvin_gnu', 'PEA_P1_M.f45_g37_rx1.A.melvin_gnu']

    >>> get_full_test_names(["acme_tiny", "^NCK.f19_g16_rx1.A"], "melvin", "gnu")
    ['ERS.f19_g16_rx1.A.melvin_gnu']
    """
    acme_test_suites = get_test_suites()

    tests_to_run = set()
    negations = set()

    for testarg in testargs:
        if (testarg.startswith("^")):
            negations.add(testarg[1:])
        elif (testarg in acme_test_suites):
            tests_to_run.update(get_test_suite(testarg, machine, compiler))
        else:
            tests_to_run.add(acme_util.get_full_test_name(testarg, machine, compiler))

    for negation in negations:
        if (negation in acme_test_suites):
            for test, testmod in get_test_suite(negation, machine, compiler):
                fullname = acme_util.get_full_test_name(test, machine, compiler, testmod)
                if (fullname in tests_to_run):
                    tests_to_run.remove(fullname)
        else:
            fullname = acme_util.get_full_test_name(negation, machine, compiler)
            if (fullname in tests_to_run):
                tests_to_run.remove(fullname)

    return list(sorted(tests_to_run))

###############################################################################
def find_all_supported_platforms():
###############################################################################
    """
    Returns a set of all ACME supported platforms as defined in the
    XML configuration file config_machines.xml in the ACME source
    tree. A platform is defined by a triple (machine name, compiler,
    mpi library).
    """
    machines = acme_util.get_machines()
    platform_set = set()

    for machine in machines:
        compilers, mpilibs = acme_util.get_machine_info(["COMPILERS", "MPILIBS"], machine=machine)
        for compiler in compilers:
            for mpilib in mpilibs:
                platform_set.add((machine, compiler, mpilib))

    return list(platform_set)

###############################################################################
def find_all_platforms(xml_file):
###############################################################################
    f = open(xml_file, "r")
    lines = f.readlines()
    f.close()
    platform_set = set()

    for line in lines:
        if "<machine" in line:
            i1 = line.index("compiler") + len('compiler="')
            i2 = line.index('"', i1)
            compiler = line[i1:i2]
            j1 = line.index(">") + 1
            j2 = line.index("<", j1)
            machine = line[j1:j2]
            platform_set.add((machine, compiler))

    return list(platform_set)

###############################################################################
def generate_acme_test_entries(category, platforms):
###############################################################################
    test_file = tempfile.NamedTemporaryFile(mode="w", delete = False)

    for machine, compiler in platforms:
        tests = get_test_suite(category, machine, compiler)
        test_file.write("\n".join(tests))

    name = test_file.name
    test_file.close()
    return name

###############################################################################
def update_acme_tests(xml_file, categories, platform=None):
###############################################################################
    # Retrieve all supported ACME platforms, killing the third entry (MPI lib)
    # for the moment.
    supported_platforms = [p[:2] for p in find_all_supported_platforms()]

    # Fish all of the existing machine/compiler combos out of the XML file.
    if (platform is not None):
        platforms = [tuple(platform.split(","))]
    else:
        platforms = find_all_platforms(xml_file)
        # Prune the non-supported platforms from our list.
        for p in platforms:
            if p not in supported_platforms:
                acme_util.verbose_print("pruning unsupported platform %s"%repr(p))
        platforms = [p for p in platforms if p in supported_platforms]

    manage_xml_entries = os.path.join(acme_util.get_cime_root(), "scripts", "manage_testlists")

    expect(os.path.isfile(manage_xml_entries),
           "Couldn't find manage_testlists, expected it to be here: '%s'" % manage_xml_entries)

    for category in categories:
        # Remove any existing acme test category from the file.
        if (platform is None):
            acme_util.run_cmd("%s -component allactive -removetests -category %s" % (manage_xml_entries, category))
        else:
            acme_util.run_cmd("%s -component allactive -removetests -category %s -machine %s -compiler %s"
                              % (manage_xml_entries, category, platforms[0][0], platforms[0][1]))

        # Generate a list of test entries corresponding to our suite at the top
        # of the file.
        new_test_file = generate_acme_test_entries(category, platforms)
        acme_util.run_cmd("%s -component allactive -addlist -file %s -category %s" %
                          (manage_xml_entries, new_test_file, category))
        os.unlink(new_test_file)

    print "SUCCESS"
