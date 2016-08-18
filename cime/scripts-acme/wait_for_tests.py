import os, time, threading, Queue, socket, signal, distutils.spawn, shutil, glob

import xml.etree.ElementTree as xmlet

import acme_util
from acme_util import expect, warning, verbose_print
from collections import OrderedDict

TEST_STATUS_FILENAME      = "TestStatus"
TEST_PENDING_STATUS       = "PEND"
TEST_PASS_STATUS          = "PASS"
TEST_FAIL_STATUS          = "FAIL"
TEST_DIFF_STATUS          = "DIFF"
NAMELIST_FAIL_STATUS      = "NLFAIL"
COMMENT_STATUS            = "COMMENT"
SLEEP_INTERVAL_SEC        = .1
THROUGHPUT_TEST_STR       = "tputcomp"
MEMORY_TEST_STR           = "memcomp"
SIGNAL_RECEIVED           = False
ACME_MAIN_CDASH           = "ACME_Climate"
CDASH_DEFAULT_BUILD_GROUP = "ACME_Latest"

NAMELIST_PHASE            = "nlcomp"
HIST_COMPARE_PHASE        = "compare"
RUN_PHASE                 = "RUN"

###############################################################################
def signal_handler(*_):
###############################################################################
    global SIGNAL_RECEIVED
    SIGNAL_RECEIVED = True

###############################################################################
def set_up_signal_handlers():
###############################################################################
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

###############################################################################
def get_test_time(test_path):
###############################################################################
    cmd = "grep TIME %s" % os.path.join(test_path, TEST_STATUS_FILENAME)
    stat, output, _ = acme_util.run_cmd(cmd, ok_to_fail=True)
    if (stat == 0):
        return int(output.split()[-1])
    else:
        warning("No timing data found in %s" % test_path)
        return 0

###############################################################################
def get_test_output(test_path):
###############################################################################
    output_file = os.path.join(test_path, "TestStatus.log")
    if (os.path.exists(output_file)):
        return open(output_file, 'r').read()
    else:
        warning("File '%s' not found" % output_file)
        return ""

###############################################################################
def create_cdash_test_xml(results, cdash_build_name, cdash_build_group, utc_time, current_time, hostname):
###############################################################################
    git_commit = acme_util.get_current_commit(repo=acme_util.get_cime_root())

    data_rel_path = os.path.join("Testing", utc_time)

    site_elem = xmlet.Element("Site")

    if ("JENKINS_START_TIME" in os.environ):
        time_info_str = "Total testing time: %d seconds" % (current_time - int(os.environ["JENKINS_START_TIME"]))
    else:
        time_info_str = ""

    site_elem.attrib["BuildName"] = cdash_build_name
    site_elem.attrib["BuildStamp"] = "%s-%s" % (utc_time, cdash_build_group)
    site_elem.attrib["Name"] = hostname
    site_elem.attrib["OSName"] = "Linux"
    site_elem.attrib["Hostname"] = hostname
    site_elem.attrib["OSVersion"] = "Commit: %s%s" % (git_commit, time_info_str)

    testing_elem = xmlet.SubElement(site_elem, "Testing")

    start_date_time_elem = xmlet.SubElement(testing_elem, "StartDateTime")
    start_date_time_elem.text = time.ctime(current_time)

    start_test_time_elem = xmlet.SubElement(testing_elem, "StartTestTime")
    start_test_time_elem.text = str(int(current_time))

    test_list_elem = xmlet.SubElement(testing_elem, "TestList")
    for test_name in sorted(results):
        test_elem = xmlet.SubElement(test_list_elem, "Test")
        test_elem.text = test_name

    for test_name in sorted(results):
        test_path, test_status = results[test_name]
        test_passed = test_status == TEST_PASS_STATUS
        test_norm_path = test_path if os.path.isdir(test_path) else os.path.dirname(test_path)

        full_test_elem = xmlet.SubElement(testing_elem, "Test")
        if (test_passed):
            full_test_elem.attrib["Status"] = "passed"
        elif (test_status == NAMELIST_FAIL_STATUS):
            full_test_elem.attrib["Status"] = "notrun"
        else:
            full_test_elem.attrib["Status"] = "failed"

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
            ("text/string",    "Completion Status", "Not Completed" if test_status == TEST_PENDING_STATUS else "Completed"),
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

###############################################################################
def create_cdash_upload_xml(results, cdash_build_name, cdash_build_group, utc_time, hostname):
###############################################################################

    data_rel_path = os.path.join("Testing", utc_time)

    try:
        log_dir = "%s_logs" % cdash_build_name

        need_to_upload = False

        for test_name, test_data in results.iteritems():
            test_path, test_status = test_data

            if (test_status not in [TEST_PASS_STATUS, NAMELIST_FAIL_STATUS]):
                full_results = parse_test_status_file(test_path)[0]

                if ("BUILD" in full_results): # If did not even make it to build phase, no useful logs
                    if ( full_results["BUILD"] != TEST_PASS_STATUS or
                         ("RUN" in full_results and full_results["RUN"] != TEST_PASS_STATUS) ):

                        param = "EXEROOT" if full_results["BUILD"] != TEST_PASS_STATUS else "RUNDIR"
                        src_dir = acme_util.run_cmd("./xmlquery %s -value" % param, from_dir=os.path.dirname(test_path))
                        log_dst_dir = os.path.join(log_dir, "%s_%s_logs" % (test_name, param))
                        os.makedirs(log_dst_dir)
                        for log_file in glob.glob(os.path.join(src_dir, "*log*")):
                            shutil.copy(log_file, log_dst_dir)

                        need_to_upload = True

        if (need_to_upload):

            tarball = "%s.tar.gz" % log_dir
            if (os.path.exists(tarball)):
                os.remove(tarball)

            acme_util.run_cmd("tar -cf - %s | gzip -c > %s" % (log_dir, tarball))
            base64 = acme_util.run_cmd("base64 %s" % tarball)

            xml_text = \
r"""<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" href="Dart/Source/Server/XSL/Build.xsl <file:///Dart/Source/Server/XSL/Build.xsl> "?>
<Site BuildName="%s" BuildStamp="%s-%s" Name="%s" Generator="ctest3.0.0">
<Upload>
<File filename="%s">
<Content encoding="base64">
%s
</Content>
</File>
</Upload>
</Site>
""" % (cdash_build_name, utc_time, cdash_build_group, hostname, os.path.abspath(tarball), base64)

            with open(os.path.join(data_rel_path, "Upload.xml"), "w") as fd:
                fd.write(xml_text)

    finally:
        if (os.path.isdir(log_dir)):
            shutil.rmtree(log_dir)

###############################################################################
def create_cdash_xml(results, cdash_build_name, cdash_project, cdash_build_group):
###############################################################################

    #
    # Create dart config file
    #

    current_time = time.time()

    utc_time_tuple = time.gmtime(current_time)
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
DropLocation: /submit.php?project=%s
DropSiteUser:
DropSitePassword:
DropSiteMode:
DropMethod: http
TriggerSite:
ScpCommand: %s

# Dashboard start time
NightlyStartTime: %s UTC
""" % (os.getcwd(), os.getcwd(), hostname,
       cdash_build_name, cdash_project, distutils.spawn.find_executable("scp"), cdash_timestamp)

    with open("DartConfiguration.tcl", "w") as dart_fd:
        dart_fd.write(dart_config)

    utc_time = time.strftime('%Y%m%d-%H%M', utc_time_tuple)
    os.makedirs(os.path.join("Testing", utc_time))

    # Make tag file
    with open("Testing/TAG", "w") as tag_fd:
        tag_fd.write("%s\n%s\n" % (utc_time, cdash_build_group))

    create_cdash_test_xml(results, cdash_build_name, cdash_build_group, utc_time, current_time, hostname)

    create_cdash_upload_xml(results, cdash_build_name, cdash_build_group, utc_time, hostname)

    acme_util.run_cmd("ctest -VV -D NightlySubmit", verbose=True)

###############################################################################
def reduce_stati(stati, check_throughput=False, check_memory=False, ignore_namelists=False):
###############################################################################
    """
    Given a collection of stati for a test, produce a single result. Preference
    is given to unfinished stati since we don't want to stop waiting for a test
    that hasn't finished. Namelist diffs are given the lowest precedence.
    """
    rv = TEST_PASS_STATUS
    for phase, status in stati.iteritems():
        if (status == TEST_PENDING_STATUS):
            return status

        elif (status != TEST_PASS_STATUS):
            if ( (not check_throughput and THROUGHPUT_TEST_STR in phase) or
                 (not check_memory and MEMORY_TEST_STR in phase) or
                 (ignore_namelists and phase == NAMELIST_PHASE) ):
                continue

            if (status == NAMELIST_FAIL_STATUS):
                if (rv == TEST_PASS_STATUS):
                    rv = NAMELIST_FAIL_STATUS

            elif (rv in [NAMELIST_FAIL_STATUS, TEST_PASS_STATUS] and phase == HIST_COMPARE_PHASE):
                rv = TEST_DIFF_STATUS

            else:
                rv = status

    return rv

###############################################################################
def parse_test_status(file_contents):
###############################################################################
    """
    Returns {ordered dict of phase->status}, test_name
    """
    rv = OrderedDict()
    test_name = None
    for line in file_contents.splitlines():
        if (line and line.split()[0] == COMMENT_STATUS):
            pass # skip comments and lines that are blank
        elif (len(line.split()) == 3):
            status, curr_test_name, phase = line.split()
            if (test_name is None):
                test_name = curr_test_name
            if (phase in rv):
                # Phase names don't matter here, just need something unique
                rv[phase] = reduce_stati({"%s_" % phase : status, phase : rv[phase]})
            else:
                rv[phase] = status
        else:
            warning("In TestStatus file for test '%s', line '%s' not in expected format" % (test_name, line))

    return rv, test_name

###############################################################################
def parse_test_status_file(file_name):
###############################################################################
    with open(file_name, "r") as fd:
        return parse_test_status(fd.read())

###############################################################################
def interpret_status(file_contents, check_throughput=False, check_memory=False, ignore_namelists=False):
###############################################################################
    r"""
    >>> interpret_status('PASS testname RUN')
    ('testname', 'PASS')
    >>> interpret_status('PASS testname BUILD\nPEND testname RUN')
    ('testname', 'PEND')
    >>> interpret_status('FAIL testname BUILD\nPEND testname RUN')
    ('testname', 'PEND')
    >>> interpret_status('PASS testname BUILD\nPASS testname RUN')
    ('testname', 'PASS')
    >>> interpret_status('PASS testname RUN\nFAIL testname tputcomp')
    ('testname', 'PASS')
    >>> interpret_status('PASS testname RUN\nFAIL testname tputcomp', check_throughput=True)
    ('testname', 'FAIL')
    >>> interpret_status('PASS testname RUN\nNLFAIL testname nlcomp')
    ('testname', 'NLFAIL')
    >>> interpret_status('PASS testname RUN\nNLFAIL testname nlcomp', ignore_namelists=True)
    ('testname', 'PASS')
    >>> interpret_status('PASS testname compare\nNLFAIL testname nlcomp\nFAIL testname compare')
    ('testname', 'DIFF')
    """
    statuses, test_name = parse_test_status(file_contents)
    reduced_status = reduce_stati(statuses, check_throughput, check_memory, ignore_namelists)

    if (RUN_PHASE not in statuses.keys() and reduced_status != TEST_FAIL_STATUS):
        warning("Very odd: Waiting for test '%s' that has no run phase but did not fail?!?!" % test_name)

    return test_name, reduced_status

###############################################################################
def interpret_status_file(file_name, check_throughput=False, check_memory=False, ignore_namelists=False):
###############################################################################
    with open(file_name, "r") as fd:
        return interpret_status(fd.read(), check_throughput, check_memory, ignore_namelists)

###############################################################################
def wait_for_test(test_path, results, wait, check_throughput, check_memory, ignore_namelists):
###############################################################################
    if (os.path.isdir(test_path)):
        test_status_filepath = os.path.join(test_path, TEST_STATUS_FILENAME)
    else:
        test_status_filepath = test_path
    verbose_print("Watching file: '%s'" % test_status_filepath)

    while (True):
        if (os.path.exists(test_status_filepath)):
            test_name, test_status = interpret_status_file(test_status_filepath, check_throughput, check_memory, ignore_namelists)

            if (test_status == TEST_PENDING_STATUS and (wait and not SIGNAL_RECEIVED)):
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
def wait_for_tests_impl(test_paths, no_wait=False, check_throughput=False, check_memory=False, ignore_namelists=False):
###############################################################################
    results = Queue.Queue()

    for test_path in test_paths:
        t = threading.Thread(target=wait_for_test, args=(test_path, results, not no_wait, check_throughput, check_memory, ignore_namelists))
        t.daemon = True
        t.start()

    while threading.active_count() > 1:
        time.sleep(1)

    test_results = {}
    completed_test_paths = []
    while (not results.empty()):
        test_name, test_path, test_status = results.get()
        if (test_name in test_results):
            prior_path, prior_status = test_results[test_name]
            if (test_status == prior_status):
                warning("Test name '%s' was found in both '%s' and '%s'" %
                        (test_name, test_path, prior_path))
            else:
                raise SystemExit("Test name '%s' was found in both '%s' and '%s' with different results" %
                                 (test_name, test_path, prior_path))

        test_results[test_name] = (test_path, test_status)
        completed_test_paths.append(test_path)

    expect(set(test_paths) == set(completed_test_paths),
           "Missing results for test paths: %s" % (set(test_paths) - set(completed_test_paths)) )

    return test_results

###############################################################################
def wait_for_tests(test_paths,
                   no_wait=False,
                   check_throughput=False,
                   check_memory=False,
                   ignore_namelists=False,
                   cdash_build_name=None,
                   cdash_project=ACME_MAIN_CDASH,
                   cdash_build_group=CDASH_DEFAULT_BUILD_GROUP):
###############################################################################
    # Set up signal handling, we want to print results before the program
    # is terminated
    set_up_signal_handlers()

    test_results = wait_for_tests_impl(test_paths, no_wait, check_throughput, check_memory, ignore_namelists)

    all_pass = True
    for test_name, test_data in sorted(test_results.iteritems()):
        test_path, test_status = test_data
        print "Test '%s' finished with status '%s'" % (test_name, test_status)
        verbose_print("    Path: %s" % test_path)
        all_pass &= test_status == TEST_PASS_STATUS

    if (cdash_build_name):
        create_cdash_xml(test_results, cdash_build_name, cdash_project, cdash_build_group)

    return all_pass
