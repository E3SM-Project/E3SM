import os, tempfile

import acme_util
from acme_util import expect, warning

# Here are the tests belonging to acme suites. Format is
# <test>.<grid>.<compset>.
# suite_name -> (inherits_from, [(test, mods)])
_TEST_SUITES = {
    "acme_tiny" : (None,
                   [("ERS.f19_g16_rx1.A", None),
                    ("NCK.f19_g16_rx1.A", None)]
                   ),
    "acme_developer" : (None,
                        [("ERS.f19_g16_rx1.A", None),
                         ("ERS.f45_g37.B1850C5", None),
                         ("ERS.f45_g37_rx1.DTEST", None),
                         ("ERS.ne30_g16_rx1.A", None),
                         ("ERS_IOP.f19_g16_rx1.A", None),
                         ("ERS_IOP.f45_g37_rx1.DTEST", None),
                         ("ERS_IOP.ne30_g16_rx1.A", None),
                         ("ERS_IOP4c.f19_g16_rx1.A", None),
                         ("ERS_IOP4c.ne30_g16_rx1.A", None),
                         ("ERS_IOP4p.f19_g16_rx1.A", None),
                         ("ERS_IOP4p.ne30_g16_rx1.A", None),
                         ("NCK.f19_g16_rx1.A", None),
                         ("PEA_P1_M.f45_g37_rx1.A", None),
                         ("SMS.ne30_f19_g16_rx1.A", None),
                         ("SMS.f19_f19.I1850CLM45CN", None),
                         ("SMS.f09_g16.I1850CLM45CN", None),
                         ("SMS.hcru_hcru.I1850CRUCLM45CN", None),
                         ("ERS_Ld5.T62_mpas120.C_MPAS_NORMAL_YEAR", None),
                         ("ERS.f09_g16_g.MPASLI_ONLY", None),
                         ("SMS.T62_mpas120_gis20.MPAS_LISIO_TEST", None),
                         ("SMS.f09_g16_a.IGCLM45_MLI", None)]
                        ),
    "acme_integration" : ("acme_developer",
                          [("ERS.ne30_g16.B1850C5", None),
                           ("ERS.f19_f19.FAMIPC5", None),
                           ("ERS.ne16_g37.B1850C5", None),
                           ("ERS_D.f45_g37.B1850C5", None),
                           ("ERS_IOP_Ld3.f19_f19.FAMIPC5", None),
                           ("ERS_Ld3.ne16_g37.FC5", None),
                           ("ERS_Ld3.ne30_ne30.FC5", None),
                           ("ERT.ne16_g37.B1850C5", None),
                           ("PET_PT.f19_g16.X", None),
                           ("PET_PT.f45_g37_rx1.A", None),
                           ("PFS.ne30_ne30.FC5", None),
                           ("SEQ_IOP_PFC.f19_g16.X", None),
                           ("SEQ_PFC.f45_g37.B1850C5", None),
                           ("SMS.ne16_ne16.FC5AQUAP", None),
                           ("SMS_D.f19_g16.B20TRC5", None),
                           ("SMS_D_Ld3.ne16_ne16.FC5", None)]
                          ),
}

###############################################################################
def get_test_suite(suite):
###############################################################################
    expect(suite in _TEST_SUITES, "Unknown test suite: '%s'" % suite)
    inherits_from, tests = _TEST_SUITES[suite]
    if (inherits_from is not None):
        inherited_tests = get_test_suite(inherits_from)
        expect(len(set(tests) & set(inherited_tests)) == 0,
               "Tests %s defined in multiple suites" % ", ".join(set(tests) & set(inherited_tests)))
        tests.extend(inherited_tests)
    return tests

###############################################################################
def get_test_suites():
###############################################################################
    return _TEST_SUITES.keys()

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

    tests = get_test_suite(category)
    for test, mods in tests:
        for machine, compiler in platforms:
            test_file.write("%s.%s_%s%s\n"%(test, machine, compiler, "" if mods is None else ".%s" % mods))

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
            acme_util.run_cmd("%s -component allactive -removetests -category %s" % (manage_xml_entries, category), verbose=True)
        else:
            acme_util.run_cmd("%s -component allactive -removetests -category %s -machine %s -compiler %s"
                              % (manage_xml_entries, category, platforms[0][0], platforms[0][1]), verbose=True)

        # Generate a list of test entries corresponding to our suite at the top
        # of the file.
        new_test_file = generate_acme_test_entries(category, platforms)
        acme_util.run_cmd("%s -component allactive -addlist -file %s -category %s" %
                          (manage_xml_entries, new_test_file, category), verbose=True)
        os.unlink(new_test_file)

    print "SUCCESS"
