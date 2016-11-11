import os, tempfile, logging
import CIME.utils
from CIME.utils import expect, run_cmd_no_fail
from CIME.XML.machines import Machines

# Here are the tests belonging to acme suites. Format is
# <test>.<grid>.<compset>.
# suite_name -> (inherits_from, [test [, mods[, machines]]])
#   To elaborate, if no mods are needed, a string representing the testname is all that is needed.
#   If testmods are needed, a 2-ple must be provided  (test, mods)
#   If you want to restrict the test mods to certain machines, than a 3-ple is needed (test, mods, [machines])
_TEST_SUITES = {
    "cime_tiny" : (None,
                   ("ERS.f19_g16_rx1.A",
                    "NCK.f19_g16_rx1.A")
                   ),

    "cime_test_only_pass" : (None,
                   ("TESTRUNPASS_Mmpi-serial.f19_g16_rx1.A",
                    "TESTRUNPASS_Mmpi-serial.ne30_g16_rx1.A",
                    "TESTRUNPASS_Mmpi-serial.f45_g37_rx1.A")
                   ),

    "cime_test_only_slow_pass" : (None,
                   ("TESTRUNSLOWPASS_Mmpi-serial.f19_g16_rx1.A",
                    "TESTRUNSLOWPASS_Mmpi-serial.ne30_g16_rx1.A",
                    "TESTRUNSLOWPASS_Mmpi-serial.f45_g37_rx1.A")
                   ),

    "cime_test_only" : (None,
                   ("TESTBUILDFAIL.f19_g16_rx1.A",
                    "TESTBUILDFAILEXC.f19_g16_rx1.A",
                    "TESTRUNFAIL_Mmpi-serial.f19_g16_rx1.A",
                    "TESTRUNFAILEXC_Mmpi-serial.f19_g16_rx1.A",
                    "TESTRUNPASS_Mmpi-serial.f19_g16_rx1.A",
                    "TESTMEMLEAKFAIL_Mmpi-serial.f19_g16.X",
                    "TESTMEMLEAKPASS_Mmpi-serial.f19_g16.X")
                   ),

    "cime_developer" : (None,
                            ("NCK_Ld3.f45_g37_rx1.A",
                             "ERI.f45_g37.X",
                             "SEQ_Ln9.f19_g16_rx1.A",
                             "ERS_Ld3.ne30_g16_rx1.A",
                             "ERS_N2_Ld3.f19_g16_rx1.A",
                             "ERR_Ld3.f45_g37_rx1.A",
                             "SMS_D_Ln9.f19_g16_rx1.A")
                            ),

    #
    # ACME tests below
    #

    "acme_runoff_developer" : (None,
                             ("ERS.f19_f19.IM1850CLM45CN",
                              "ERS.f19_f19.IMCLM45")
                             ),

    "acme_land_developer" : ("acme_runoff_developer",
                             ("ERS.f19_f19.I1850CLM45CN",
                              "ERS.f09_g16.I1850CLM45CN",
                              "SMS.hcru_hcru.I1850CRUCLM45CN",
                             ("SMS_Ly3.1x1_smallvilleIA.ICLM45CNCROP", "force_netcdf_pio"),
                              "ERS.ne11_oQU240.I20TRCLM45",
                              "ERS.f09_g16.IMCLM45BC")
                             ),

    "acme_atm_developer" : (None,
                            ("ERS.ne16_ne16.FC5MAM4",
                             "ERS.ne16_ne16.FC5PLMOD",
                             "ERS.ne16_ne16.FC5CLBMG2",
                             "ERS.ne16_ne16.FC5CLBMG2MAM4",
                             "ERS.ne16_ne16.FC5CLBMG2MAM4MOM",
                             "ERS.ne16_ne16.FC5CLBMG2MAM4RESUS",
                             "ERS.ne16_ne16.FC5CLBMG2LINMAM4RESUSMOM",
                             "ERS.f19_g16.FC5CLBMG2MAM4RESUSBC",
                             "ERS.f19_g16.FC5CLBMG2MAM4RESUSMOMBC",
                             "ERS.f19_g16.FC5ATMMOD",
                             "ERS_Ld5.ne16_ne16.FC5ATMMODCOSP",
                             "SMS.f19_g16.FC5ATMMOD",
                             "SMS.f19_g16.FC5ATMMODCOSP",
                             "SMS_D.f19_g16.FC5ATMMODCOSP")
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
                         "HOMME_P24.f19_g16_rx1.A",
                         "NCK.f19_g16_rx1.A",
                         "SMS.ne30_f19_g16_rx1.A",
                         "ERS_Ld5.T62_oQU120.CMPASO-NYF",
                         "ERS.f09_g16_g.MPASLI_ONLY",
                         "ERS_Ld5.ne16_ne16.FC5ATMMODCOSP",
                         "SMS.T62_oQU120_ais20.MPAS_LISIO_TEST",
                         "SMS.f09_g16_a.IGCLM45_MLI",
                         "SMS_D_Ln5.ne16_ne16.FC5AV1F",
                         "SMS_D_Ln5.ne16_ne16.FC5AV1C",
                         "SMS_D_Ln5.ne16_ne16.FC5AV1C-01",
                         "SMS_D_Ln5.ne16_ne16.FC5AV1C-02",
                         "SMS_D_Ln5.ne16_ne16.FC5AV1C-03",
                         "SMS_D_Ln5.ne16_ne16.FC5AV1C-04",
                         "SMS_D_Ln1.ne30_ne30.FC5AV1C-04",
                         "SMS_D_Ln1.ne30_oEC.F1850C5AV1C-02",
                         "SMS_D_Ln5.ne16_ne16.F1850C5AV1C-04",
                         "SMS_D_Ln5.ne16_ne16.F20TRC5AV1C-03",
                         "SMS_D_Ln5.ne16_ne16.FC5AV1C-04P",
                         "SMS_D_Ld1.ne16_ne16.FC5ATMMOD")
                        ),

    "acme_integration" : ("acme_developer",
                          ("ERS.ne11_oQU240.A_WCYCL1850",
                           "ERS.f19_f19.FAMIPC5",
                           "ERS.ne16_ne16.FC5PM",
                           "ERS.ne16_ne16.FC5PLMOD",
                           "ERS.ne16_ne16.FC5MAM4",
                           "ERS_IOP_Ld3.f19_f19.FAMIPC5",
                           "ERS_Ld3.ne16_g37.FC5",
                           "ERS_Ld3.ne30_ne30.FC5",
                          #"ERT_Ld31.ne16_g37.B1850C5",#add this line back in with the new correct compset
                           ("PET_Ln9.ne30_ne30.FC5", "cam-outfrq9s"),
                           "PET.f19_g16.X",
                           "PET.f45_g37_rx1.A",
                           "PET_Ln9.ne30_oEC.A_WCYCL2000",
                           "ERP_Ln3.ne30_oEC.A_WCYCL2000",
                           "PFS.ne30_ne30.FC5",
                           "SEQ_IOP.f19_g16.X",
                           "SMS.ne30_oEC.A_WCYCL2000",
                           "SMS.ne16_ne16.FC5AQUAP",
                           "SMS_D_Ld3.ne16_ne16.FC5",
                           "SMS.f09_g16_a.MPASLIALB_ONLY",
                           "ERS.ne16_ne16.FC5ATMMOD",
                           "ERS_Ld5.ne16_ne16.FC5AV1F",
                           "ERS_Ld5.ne16_ne16.FC5AV1C",
                           "ERS_Ld5.ne16_ne16.FC5AV1C-01",
                           "ERS_Ld5.ne16_ne16.FC5AV1C-02",
                           "ERS_Ld5.ne16_ne16.FC5AV1C-03",
                           "ERS_Ld5.ne16_ne16.FC5AV1C-04",
                           "ERS_Ld5.ne30_oEC.F1850C5AV1C-02",
                           "ERS_Ld5.ne16_ne16.F1850C5AV1C-04",
                           "ERS_Ld5.ne16_ne16.F20TRC5AV1C-03",
                           "SMS_D_Ld1.ne16_ne16.FC5ATMMODCOSP")
                          ),
}

###############################################################################
def get_test_suite(suite, machine=None, compiler=None):
###############################################################################
    """
    Return a list of FULL test names for a suite.
    """
    expect(suite in _TEST_SUITES, "Unknown test suite: '%s'" % suite)
    machobj = Machines(machine=machine)
    machine = machobj.get_machine_name()

    if(compiler is None):
        compiler = machobj.get_default_compiler()
    expect(machobj.is_valid_compiler(compiler),"Compiler %s not valid for machine %s" %
           (compiler,machine))

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

        tests.append(CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler, testmod=test_mod))

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
    acme_test_suites = get_test_suites()

    tests_to_run = set()
    negations = set()

    for testarg in testargs:
        # remove any whitespace in name
        testarg = testarg.strip()
        if (testarg.startswith("^")):
            negations.add(testarg[1:])
        elif (testarg in acme_test_suites):
            tests_to_run.update(get_test_suite(testarg, machine, compiler))
        else:
            tests_to_run.add(CIME.utils.get_full_test_name(testarg, machine=machine, compiler=compiler))

    for negation in negations:
        if (negation in acme_test_suites):
            for test, testmod in get_test_suite(negation, machine, compiler):
                fullname = CIME.utils.get_full_test_name(test, machine=machine, compiler=compiler, testmod=testmod)
                if (fullname in tests_to_run):
                    tests_to_run.remove(fullname)
        else:
            fullname = CIME.utils.get_full_test_name(negation, machine=machine, compiler=compiler)
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
    # TODO - Fix
    pass
    # machines = CIME.utils.get_machines()
    # machobj = Machines(machine=machines)
    # platform_set = set()

    # for machine in machines:
    #     machobj.set_machine(machine)
    #     compilers, mpilibs = machobj.get_value("COMPILERS"), machobj.get_value("MPILIBS")
    #     for compiler in compilers:
    #         for mpilib in mpilibs:
    #             platform_set.add((machine, compiler, mpilib))

    # return list(platform_set)

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
                logging.info("pruning unsupported platform %s"%repr(p))
        platforms = [p for p in platforms if p in supported_platforms]

    manage_xml_entries = os.path.join(CIME.utils.get_cime_root(), "scripts", "manage_testlists")

    expect(os.path.isfile(manage_xml_entries),
           "Couldn't find manage_testlists, expected it to be here: '%s'" % manage_xml_entries)

    for category in categories:
        # Remove any existing acme test category from the file.
        if (platform is None):
            run_cmd_no_fail("%s -model acme -component allactive -removetests -category %s" % (manage_xml_entries, category))
        else:
            run_cmd_no_fail("%s -model acme -component allactive -removetests -category %s -machine %s -compiler %s"
                            % (manage_xml_entries, category, platforms[0][0], platforms[0][1]))

        # Generate a list of test entries corresponding to our suite at the top
        # of the file.
        new_test_file = generate_acme_test_entries(category, platforms)
        run_cmd_no_fail("%s -model acme -component allactive -addlist -file %s -category %s" %
                        (manage_xml_entries, new_test_file, category))
        os.unlink(new_test_file)

    print "SUCCESS"
