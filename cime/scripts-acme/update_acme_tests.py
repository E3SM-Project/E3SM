import os, shutil, tempfile

import acme_util
from acme_util import expect, warning

# Here are the tests belonging to acme suites. Format is
# <test>.<grid>.<compset>.
# suite_name -> (inherits_from, [testlist])
__TEST_SUITES = {
    "acme_tiny" : (None,
                   ["ERS.f19_g16_rx1.A",
                    "NCK.f19_g16_rx1.A"]
                   ),
    "acme_developer" : (None,
                        ['ERS.f19_g16_rx1.A',
                         'ERS.f45_g37.B1850C5',
                         'ERS.f45_g37_rx1.DTEST',
                         'ERS.ne30_g16_rx1.A',
                         'ERS_IOP.f19_g16_rx1.A',
                         'ERS_IOP.f45_g37_rx1.DTEST',
                         'ERS_IOP.ne30_g16_rx1.A',
                         'ERS_IOP4c.f19_g16_rx1.A',
                         'ERS_IOP4c.ne30_g16_rx1.A',
                         'ERS_IOP4p.f19_g16_rx1.A',
                         'ERS_IOP4p.ne30_g16_rx1.A',
                         'NCK.f19_g16_rx1.A',
                         'PEA_P1_M.f45_g37_rx1.A',
                         'SMS.ne30_f19_g16_rx1.A',
                         'SMS.f19_f19.I1850CLM45CN',
			 'SMS.f09_g16.I1850CLM45CN',
			 'SMS.hcru_hcru.I1850CRUCLM45CN']
                        ),
    "acme_integration" : ("acme_developer",
                          ["ERS.ne30_g16.B1850C5",
                           "ERS.f19_f19.FAMIPC5",
                           "ERS.ne16_g37.B1850C5",
                           "ERS_D.f45_g37.B1850C5",
                           "ERS_IOP_Ld3.f19_f19.FAMIPC5",
                           "ERS_Ld3.ne16_g37.FC5",
                           "ERS_Ld3.ne30_ne30.FC5",
                           "ERT.ne16_g37.B1850C5",
                           "PET_PT.f19_g16.X",
                           "PET_PT.f45_g37_rx1.A",
                           "PFS.ne30_ne30.FC5",
                           "SEQ_IOP_PFC.f19_g16.X",
                           "SEQ_PFC.f45_g37.B1850C5",
                           "SMS.ne16_ne16.FC5AQUAP",
                           "SMS_D.f19_g16.B20TRC5",
                           "SMS_D_Ld3.ne16_ne16.FC5"]
                          ),
}

###############################################################################
def get_test_suite(suite):
###############################################################################
    expect(suite in __TEST_SUITES, "Unknown test suite: '%s'" % suite)
    inherits_from, tests = __TEST_SUITES[suite]
    if (inherits_from is not None):
        inherited_tests = get_test_suite(inherits_from)
        expect(len(set(tests) & set(inherited_tests)) == 0,
               "Tests %s defined in multiple suites" % ", ".join(set(tests) & set(inherited_tests)))
        tests.extend(inherited_tests)
    return tests

###############################################################################
def get_test_suites():
###############################################################################
    return __TEST_SUITES.keys()

###############################################################################
def find_all_supported_platforms():
###############################################################################
    """
    Returns a set of all ACME supported platforms as defined in the
    XML configuration file config_machines.xml in the ACME source
    tree. A platform is defined by a triple (machine name, compiler,
    mpi library).
    """
    import xml.etree.ElementTree as ET
    config_machines_xml = os.path.join(acme_util.get_cime_root(), 'machines-acme', 'config_machines.xml')
    tree = ET.parse(config_machines_xml)
    root = tree.getroot()
    expect(root.tag == 'config_machines',
           'The given XML file is not a valid list of machine configurations.')
    platform_set = set()

    # Each child of this root is a machine entry.
    for machine in root:
        if (machine.tag == 'machine'):
            expect('MACH' in machine.attrib, 'Invalid machine entry found')
            mach_name = machine.attrib['MACH']
            expect('COMPILERS' in [item.tag for item in machine],
                   'COMPILERS entry not found in machine %s'%mach_name)
            compilers_string = machine.find('COMPILERS').text
            compilers = [compiler.strip() for compiler in compilers_string.split(',')]
            mpilibs_string = machine.find('MPILIBS').text
            mpilibs = [mpilib.strip() for mpilib in mpilibs_string.split(',')]
            for compiler in compilers:
                for mpilib in mpilibs:
                    platform_set.add((mach_name, compiler, mpilib))
        else:
            warning("Ignoring unrecognized tag: '%s'" % machine.tag)

    return list(platform_set)

###############################################################################
def find_all_platforms(xml_file):
###############################################################################
    f = open(xml_file, 'r')
    lines = f.readlines()
    f.close()
    platform_set = set()

    for line in lines:
        if '<machine' in line:
            i1 = line.index('compiler') + len('compiler="')
            i2 = line.index('"', i1)
            compiler = line[i1:i2]
            j1 = line.index('>') + 1
            j2 = line.index('<', j1)
            machine = line[j1:j2]
            platform_set.add((machine, compiler))

    return list(platform_set)

###############################################################################
def generate_acme_test_entries(category, platforms):
###############################################################################
    tests = get_test_suite(category)
    test_file = tempfile.NamedTemporaryFile(mode='w', delete = False)
    for test in tests:
        for machine, compiler in platforms:
            test_file.write('%s.%s_%s\n'%(test, machine, compiler))
    name = test_file.name
    test_file.close()
    return name

###############################################################################
def update_acme_test(xml_file, categories, platform):
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
            acme_util.verbose_print('pruning unsupported platform %s'%repr(p))
    platforms = [p for p in platforms if p in supported_platforms]

    manage_xml_entries = os.path.join(acme_util.get_cime_root(), "scripts", "manage_testlists")

    expect(os.path.isfile(manage_xml_entries),
           "Couldn't find manage_testlists, expected it to be here: '%s'" % manage_xml_entries)

    for category in categories:
        # Remove any existing acme test category from the file.
        if (platform is None):
            output = acme_util.run_cmd('%s -component allactive -removetests -category %s' % (manage_xml_entries, category), verbose=True)
        else:
            output = acme_util.run_cmd('%s -component allactive -removetests -category %s -machine %s -compiler %s'
                                       % (manage_xml_entries, category, platforms[0][0], platforms[0][1]), verbose=True)

        # Generate a list of test entries corresponding to our suite at the top
        # of the file.
        new_test_file = generate_acme_test_entries(category, platforms)
        output = acme_util.run_cmd("%s -component allactive -addlist -file %s -category %s" %
                                   (manage_xml_entries, new_test_file, category), verbose=True)
        os.unlink(new_test_file)

    print "SUCCESS"
