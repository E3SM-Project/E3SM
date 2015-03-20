#!/usr/bin/env python

import os, doctest, time, threading, Queue, socket, signal, distutils.spawn

import xml.etree.ElementTree as xmlet

import acme_util
from acme_util import expect, warning, verbose_print

TEST_STATUS_FILENAME     = "TestStatus"
TEST_NOT_FINISHED_STATUS = ["GEN", "BUILD", "RUN", "PEND"]
TEST_PASSED_STATUS       = "PASS"
SLEEP_INTERVAL_SEC       = 1
THROUGHPUT_TEST_STR      = ".tputcomp."
SIGNAL_RECEIVED          = False

###############################################################################
def signal_handler(*_):
###############################################################################
    global SIGNAL_RECEIVED
    SIGNAL_RECEIVED = True

###############################################################################
def get_test_time(test_path):
###############################################################################
    cmd = "grep 'TOT Run Time' /dev/null $(find %s -name 'ccsm_timing*') || true" % test_path
    output = acme_util.run_cmd(cmd)

    tot_time = 0.0
    for line in output.splitlines():
        if (line != "" and not line.isspace()):
            tokens = line.split()

            if (len(tokens) < 5 or tokens[1:4] != ["TOT", "Run", "Time:"]):
                warning("Line '%s' not in expected format")
                continue

            try:
                cur_time = float(tokens[4])
                tot_time += cur_time
            except ValueError:
                warning("Line '%s' not in expected format, '%s' not a valid float" % (line, tokens[4]))

    if (tot_time == 0.0):
        warning("No timing data found in %s" % test_path)

    return tot_time

###############################################################################
def get_test_output(test_path):
###############################################################################
    output_file = os.path.join(test_path, "TestStatus.out")
    if (os.path.exists(output_file)):
        return open(output_file, 'r').read()
    else:
        warning("File '%s' not found" % output_file)
        return ""

###############################################################################
def create_cdash_xml(start_time, results, cdash_build_name):
###############################################################################

    # Create dart config file
    utc_time_tuple = time.gmtime(start_time)
    cdash_timestamp = time.strftime("%H:%M:%S", utc_time_tuple)
    hostname = acme_util.probe_machine_name()
    if (hostname is None):
        hostname = socket.gethostname().split(".")[0]
        warning("Could not convert hostname '%s' into an ACME machine name" % (hostname))
    dart_config = \
"""
SourceDirectory: %s
BuildDirectory: %s

# Site is something like machine.domain, i.e. pragmatic.crd
Site: %s

# Build name is osname-revision-compiler, i.e. Linux-2.4.2-2smp-c++
BuildName: %s

# Submission information
IsCDash: TRUE
CDashVersion:
QueryCDashVersion:
DropSite: my.cdash.org
DropLocation: /submit.php?project=ACME_Climate
DropSiteUser:
DropSitePassword:
DropSiteMode:
DropMethod: http
TriggerSite:
ScpCommand: %s

# Dashboard start time
NightlyStartTime: %s UTC
""" % (os.getcwd(), os.getcwd(), hostname,
       cdash_build_name, distutils.spawn.find_executable("scp"), cdash_timestamp)

    dart_fd = open("DartConfiguration.tcl", "w")
    dart_fd.write(dart_config)
    dart_fd.close()

    # Make necessary dirs
    subdir_name = time.strftime('%Y%m%d-%H%M', utc_time_tuple)
    data_rel_path = os.path.join("Testing", subdir_name)
    os.makedirs(data_rel_path)

    # Make tag file
    tag_fd = open("Testing/TAG", "w")
    tag_fd.write("%s\nNightly" % subdir_name)
    tag_fd.close()

    #
    # Make XML
    #

    site_elem = xmlet.Element("Site")

    # It's OK to lie for most of this stuff
    site_elem.attrib["BuildName"] = cdash_build_name
    site_elem.attrib["BuildStamp"] = "%s-Nightly" % subdir_name
    site_elem.attrib["Name"] = hostname
    site_elem.attrib["Generator"] = "ctest-2.8.11.1"
    site_elem.attrib["CompilerName"] = ""
    site_elem.attrib["OSName"] = "Linux"
    site_elem.attrib["Hostname"] = hostname
    site_elem.attrib["OSRelease"] = "2.6.32-504.el6.x86_64"
    site_elem.attrib["OSVersion"] = "#1 SMP Tue Sep 16 01:56:35 EDT 2014"
    site_elem.attrib["OSPlatform"] = "x86_64"
    site_elem.attrib["Is64Bits"] = "1"
    site_elem.attrib["VendorString"] = "GenuineIntel"
    site_elem.attrib["VendorID"] = "Intel Corporation"
    site_elem.attrib["FamilyID"] = "6"
    site_elem.attrib["ModelID"] = "44"
    site_elem.attrib["ProcessorCacheSize"] = "12288"
    site_elem.attrib["NumberOfLogicalCPU"] = "16"
    site_elem.attrib["NumberOfPhysicalCPU"] = "8"
    site_elem.attrib["TotalVirtualMemory"] = "26207"
    site_elem.attrib["TotalPhysicalMemory"] = "24016"
    site_elem.attrib["LogicalProcessorsPerPhysical"] = "2"
    site_elem.attrib["ProcessorClockFrequency"] = "2394.04"

    testing_elem = xmlet.SubElement(site_elem, "Testing")

    start_date_time_elem = xmlet.SubElement(testing_elem, "StartDateTime")
    start_date_time_elem.text = time.ctime(start_time)

    start_test_time_elem = xmlet.SubElement(testing_elem, "StartTestTime")
    start_test_time_elem.text = str(int(start_time))

    test_list_elem = xmlet.SubElement(testing_elem, "TestList")
    for test_name in sorted(results):
        test_elem = xmlet.SubElement(test_list_elem, "Test")
        test_elem.text = test_name

    for test_name in sorted(results):
        test_path, test_status = results[test_name]
        test_passed = test_status == TEST_PASSED_STATUS
        test_norm_path = test_path if os.path.isdir(test_path) else os.path.dirname(test_path)

        full_test_elem = xmlet.SubElement(testing_elem, "Test")
        full_test_elem.attrib["Status"] = "passed" if test_passed else "failed"

        name_elem = xmlet.SubElement(full_test_elem, "Name")
        name_elem.text = test_name

        path_elem = xmlet.SubElement(full_test_elem, "Path")
        path_elem.text = test_norm_path

        full_name_elem = xmlet.SubElement(full_test_elem, "FullName")
        full_name_elem.text = test_name

        full_command_line_elem = xmlet.SubElement(full_test_elem, "FullCommandLine")
        # text ?

        results_elem = xmlet.SubElement(full_test_elem, "Results")

        named_measurements = (
            ("text/string",    "Exit Code",         test_status),
            ("text/string",    "Exit Value",        "0" if test_passed else "1"),
            ("numeric_double", "Execution Time",    str(get_test_time(test_norm_path))),
            ("text/string",    "Completion Status", "Not Completed" if test_status in TEST_NOT_FINISHED_STATUS else "Completed"),
            ("text/string",    "Command line",      "create_test")
        )

        for type_attr, name_attr, value in named_measurements:
            named_measurement_elem = xmlet.SubElement(results_elem, "NamedMeasurement")
            named_measurement_elem.attrib["type"] = type_attr
            named_measurement_elem.attrib["name"] = name_attr

            value_elem = xmlet.SubElement(named_measurement_elem, "Value")
            value_elem.text = value

        measurement_elem = xmlet.SubElement(results_elem, "Measurement")

        value_elem = xmlet.SubElement(measurement_elem, "Value")
        value_elem.text = get_test_output(test_norm_path)

    elapsed_time_elem = xmlet.SubElement(testing_elem, "ElapsedMinutes")
    elapsed_time_elem.text = "0" # Skip for now

    etree = xmlet.ElementTree(site_elem)

    etree.write(os.path.join(data_rel_path, "Test.xml"))

    acme_util.run_cmd("ctest -D NightlySubmit", verbose=True)

###############################################################################
def parse_test_status_file(file_contents, status_file_path, check_throughput):
###############################################################################
    r"""
    >>> parse_test_status_file('PASS testname', '', False)
    ('testname', 'PASS')
    >>> parse_test_status_file('PASS testname \nGEN testname2', '', False)
    ('testname', 'GEN')
    >>> parse_test_status_file('PASS testname\nPASS testname2', '', False)
    ('testname', 'PASS')
    >>> parse_test_status_file('PASS testname\nFAIL testname2.tputcomp.foo', '', False)
    ('testname', 'PASS')
    >>> parse_test_status_file('PASS testname\nFAIL testname2.tputcomp.foo', '', True)
    ('testname', 'FAIL')
    """
    real_test_name = None
    for line in file_contents.splitlines():
        if (len(line.split()) == 2):
            status, test_name = line.split()
            if (real_test_name is None):
                real_test_name = test_name # just take the first one

            verbose_print("Test: '%s' has status '%s'" % (test_name, status))

            # A non-pass is OK if the failure is due to throughput and we
            # aren't checking throughput
            if (status != TEST_PASSED_STATUS and not
                (not check_throughput and THROUGHPUT_TEST_STR in test_name)):
                return real_test_name, status
        else:
            warning("In '%s', line '%s' not in expected format" % (status_file_path, line))

    if (real_test_name is None):
        warning("Empty status file: %s" % status_file_path)

    return real_test_name, TEST_PASSED_STATUS

###############################################################################
def wait_for_test(test_path, results, wait, check_throughput):
###############################################################################
    if (os.path.isdir(test_path)):
        test_status_filepath = os.path.join(test_path, TEST_STATUS_FILENAME)
    else:
        test_status_filepath = test_path
    verbose_print("Watching file: '%s'" % test_status_filepath)

    while (True):
        if (os.path.exists(test_status_filepath)):
            test_status_fd = open(test_status_filepath, "r")
            test_status_contents = test_status_fd.read()
            test_name, test_status = parse_test_status_file(test_status_contents, test_status_filepath, check_throughput)

            if (test_status in TEST_NOT_FINISHED_STATUS and (wait and not SIGNAL_RECEIVED)):
                time.sleep(SLEEP_INTERVAL_SEC)
                verbose_print("Waiting for test to finish")
            else:
                results.put( (test_name, test_path, test_status) )
                break

        else:
            if (wait and not SIGNAL_RECEIVED):
                verbose_print("File '%s' does not yet exist" % test_status_filepath)
                time.sleep(SLEEP_INTERVAL_SEC)
            else:
                test_name = os.path.abspath(test_status_filepath).split("/")[-2]
                results.put( (test_name, test_path, "File '%s' doesn't exist" % test_status_filepath) )
                break

###############################################################################
def wait_for_tests(test_paths, no_wait, check_throughput, cdash_build_name):
###############################################################################
    # Set up signal handling, we want to print results before the program
    # is terminated
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    if (cdash_build_name):
        start_time = time.time()

    results = Queue.Queue()

    for test_path in test_paths:
        t = threading.Thread(target=wait_for_test, args=(test_path, results, not no_wait, check_throughput))
        t.daemon = True
        t.start()

    while threading.active_count() > 1:
        time.sleep(SLEEP_INTERVAL_SEC)

    tests_with_results = dict()
    completed_test_paths = []
    while (not results.empty()):
        test_name, test_path, test_status = results.get()
        if (test_name in tests_with_results):
            prior_path, prior_status = tests_with_results[test_name]
            if (test_status == prior_status):
                warning("Test name '%s' was found in both '%s' and '%s'" %
                        (test_name, test_path, prior_path))
            else:
                raise SystemExit("Test name '%s' was found in both '%s' and '%s' with different results" %
                                 (test_name, test_path, prior_path))

        tests_with_results[test_name] = (test_path, test_status)
        completed_test_paths.append(test_path)

    expect(set(test_paths) == set(completed_test_paths),
           "Missing results for test paths: %s" % (set(test_paths) - set(completed_test_paths)) )

    all_pass = True
    for test_name, test_data in sorted(tests_with_results.iteritems()):
        test_path, test_status = test_data
        print "Test '%s' finished with status '%s'" % (test_name, test_status)
        print "    Path: %s" % test_path
        all_pass &= test_status == TEST_PASSED_STATUS

    if (cdash_build_name):
        create_cdash_xml(start_time, tests_with_results, cdash_build_name)

    return all_pass

###############################################################################
def run_unit_tests():
###############################################################################
    doctest.testmod()
