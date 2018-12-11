import CIME.utils
from CIME.utils import expect, convert_to_seconds
from CIME.XML.machines import Machines

# Here are the tests belonging to acme suites. Format is
# <test>.<grid>.<compset>.
# suite_name -> (inherits_from, timelimit, [test [, mods[, machines]]])
#   To elaborate, if no mods are needed, a string representing the testname is all that is needed.
#   If testmods are needed, a 2-ple must be provided  (test, mods)
#   If you want to restrict the test mods to certain machines, than a 3-ple is needed (test, mods, [machines])
_TEST_SUITES = {
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
                    "TESTRUNFAILEXC_P1.f19_g16_rx1.A",
                    "TESTRUNPASS_P1.f19_g16_rx1.A",
                    "TESTTESTDIFF_P1.f19_g16_rx1.A",
                    "TESTMEMLEAKFAIL_P1.f19_g16.X",
                    "TESTMEMLEAKPASS_P1.f19_g16.X")
                   ),

    "cime_developer" : (None, "0:10:00",
                            ("NCK_Ld3.f45_g37_rx1.A",
                             "ERI.f45_g37.X",
                             "ERIO.f45_g37.X",
                             "SEQ_Ln9.f19_g16_rx1.A",
                             "ERS.ne30_g16_rx1.A",
                             "ERS_N2.f19_g16_rx1.A",
                             "ERR.f45_g37_rx1.A",
                             "ERP.f45_g37_rx1.A",
                             "SMS_D_Ln9.f19_g16_rx1.A",
                             "DAE.f19_f19.A",
                             "PET_P32.f19_f19.A",
                             "SMS.T42_T42.S",
                             "PRE.f45_g37_rx1.ADESP")
                            ),

    #
    # ACME tests below
    #

    "acme_runoff_developer" : (None, "0:45:00",
                             ("ERS.f19_f19.IM1850CLM45CN",
                              "ERS.f19_f19.IMCLM45")
                             ),

    "acme_land_developer" : ("acme_runoff_developer", "0:45:00",
                             ("ERS.f19_f19.I1850CLM45CN",
                              "ERS.f09_g16.I1850CLM45CN",
                              "ERS.f19_f19.I20TRCLM45CN",
                              "SMS_Ld1.hcru_hcru.I1850CRUCLM45CN",
                             ("ERS.f19_g16.I1850CNECACNTBC" ,"clm-eca"),
                             ("ERS.f19_g16.I1850CNECACTCBC" ,"clm-eca"),
                             ("SMS_Ly2_P1x1.1x1_smallvilleIA.ICLM45CNCROP", "force_netcdf_pio"),
                             ("SMS_Ld4.f45_f45.ICLM45ED","clm-fates"),
                             ("ERP_Ld3.f45_f45.ICLM45ED","clm-fates"),
                             ("ERS.f19_g16.I1850CLM45","clm-betr"),
                             ("ERS.f19_g16.I1850CLM45","clm-vst"),
                             ("ERS.f09_g16.I1850CLM45CN","clm-bgcinterface"),
                              "ERS.ne11_oQU240.I20TRCLM45",
                              "ERS.f09_g16.IMCLM45BC")
                             ),

    "acme_atm_developer" : (None, None,
                            ("SMS_D_Ln5.ne4_ne4.FC5",
                             "ERS_Ln9.ne4_ne4.FC5AV1C-L",
                             ("SMS_Ln9.ne4_ne4.FC5AV1C-L", "cam-outfrq9s"),
                             "SMS_R_Ld5.T42_T42.FSCM5A97",
                             "SMS_D_Ln5.ne4_ne4.FC5AV1C-L")
                            ),

    "acme_atm_integration" : (None, None,
                              ("ERS_Ln9.ne4_ne4.FC5" ,
                               "ERP_Ln9.ne4_ne4.FC5AQUAP",
                               "PET_Ln5.ne4_ne4.FC5AV1C-L",
                               "PEM_Ln5.ne4_ne4.FC5AV1C-L",
                               ("SMS_D_Ln5.ne4_ne4.FC5AV1C-L", "cam-cosplite_nhtfrq5"),
                               "REP_Ln5.ne4_ne4.FC5AV1C-L")
                              ),

    "acme_atm_weekly" : (None, None,
                         ("SMS_Lm1.ne4_ne4.FC5AV1C-L",
                          "ERS_Ld31.ne4_ne4.FC5AV1C-L",
                          "ERP_Lm3.ne4_ne4.FC5AV1C-L",
                          "SMS_Ly1.ne4_ne4.FC5AV1C-L")
                         ),

    "acme_developer" : (("acme_land_developer","acme_atm_developer"), "0:45:00",
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
                         "HOMME_P24.f19_g16_rx1.A",
                         "NCK.f19_g16_rx1.A",
                         "SMS.ne30_f19_g16_rx1.A",
                         "ERS_Ld5.T62_oQU120.CMPASO-NYF",
                         "ERS.f09_g16_g.MPASLISIA",
                         "SMS.T62_oQU120_ais20.MPAS_LISIO_TEST",
                         "SMS.f09_g16_a.IGCLM45_MLI")
                        ),

    "acme_integration" : (("acme_developer", "acme_atm_integration"),"03:00:00",
                          ("ERS.ne11_oQU240.A_WCYCL1850",
                          #"ERT_Ld31.ne16_g37.B1850C5",#add this line back in with the new correct compset
                           "PET.f19_g16.X",
                           "PET.f45_g37_rx1.A",
                           "PET_Ln9.ne30_oECv3_ICG.A_WCYCL1850S",
                           "PEM_Ln9.ne30_oECv3_ICG.A_WCYCL1850S",
                           "ERP_Ld3.ne30_oECv3_ICG.A_WCYCL1850S",
                           "SEQ_IOP.f19_g16.X",
                           "SMS.ne30_oECv3_ICG.A_WCYCL1850S",
                           "SMS.f09_g16_a.MPASLIALB")
                          ),

    "acme_weekly" : (("acme_atm_weekly",),None,
                     () ),

    "fates" : (None, None,
                         ("ERS_Ld9.1x1_brazil.ICLM45ED",
                          "ERS_D_Ld9.1x1_brazil.ICLM45ED",
                          "SMS_D_Lm6.1x1_brazil.ICLM45ED")
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

    inherits_from, _, tests_raw = _TEST_SUITES[suite]
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
        inherits_from = [inherits_from] if isinstance(inherits_from, str) else inherits_from
        for inherits in inherits_from:
            inherited_tests = get_test_suite(inherits, machine, compiler)

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

    >>> get_recommended_test_time("PET_Ln9.ne30_ne30.FC5.skybridge_intel.cam-outfrq9s")
    '03:00:00'

    >>> get_recommended_test_time("PET_Ln20.ne30_ne30.FC5.skybridge_intel.cam-outfrq9s")
    >>>
    """
    _, _, _, _, machine, compiler, _ = CIME.utils.parse_test_name(test_full_name)
    expect(machine is not None, "%s is not a full test name" % test_full_name)

    best_time = None
    suites = get_test_suites()
    for suite in suites:
        _, rec_time, tests_raw = _TEST_SUITES[suite]
        for item in tests_raw:
            test_mod = None
            if (isinstance(item, str)):
                test_name = item
            else:
                test_name = item[0]
                if (len(item) == 2):
                    test_mod = item[1]
                else:
                    test_mod_machines = [item[2]] if isinstance(item[2], str) else item[2]
                    if (machine in test_mod_machines):
                        test_mod = item[1]

            full_test = CIME.utils.get_full_test_name(test_name, machine=machine, compiler=compiler, testmod=test_mod)

            if full_test == test_full_name and rec_time is not None:
                if best_time is None or \
                        convert_to_seconds(rec_time) < convert_to_seconds(best_time):
                    best_time = rec_time

    return best_time
