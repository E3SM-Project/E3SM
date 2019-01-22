#!/usr/bin/env python

"""
Simple script to populate CESM test database with test results.
"""

from standard_script_setup import *

from CIME.XML.env_build             import EnvBuild
from CIME.XML.env_case              import EnvCase
from CIME.XML.env_test              import EnvTest
from CIME.XML.test_reporter         import TestReporter
from CIME.utils                     import expect
from CIME.XML.generic_xml import GenericXML

import glob

###############################################################################
def parse_command_line(args):
###############################################################################
    parser = argparse.ArgumentParser()

    CIME.utils.setup_standard_logging_options(parser)

    # Parse command line options

    #parser = argparse.ArgumentParser(description='Arguements for testreporter')
    parser.add_argument("--tagname",
                        help="Name of the tag being tested.")
    parser.add_argument("--testid",
                        help="Test id, ie c2_0_a6g_ing,c2_0_b6g_gnu.")
    parser.add_argument("--testroot",
                        help="Root directory for tests to populate the database.")
    parser.add_argument("--testtype",
                        help="Type of test, prealpha or prebeta.")
    parser.add_argument("--dryrun",action="store_true",
                        help="Do a dry run, database will not be populated.")
    parser.add_argument("--dumpxml",action="store_true",
                        help="Dump XML test results to sceen.")
    args = parser.parse_args()
    CIME.utils.parse_args_and_handle_standard_logging_options(args)

    return args.testroot, args.testid, args.tagname, args.testtype, args.dryrun, args.dumpxml

###############################################################################
def get_testreporter_xml(testroot, testid, tagname, testtype):
###############################################################################
    os.chdir(testroot)

    #
    # Retrieve compiler name and mpi library
    #
    xml_file=glob.glob("*"+testid+"/env_build.xml")
    expect(len(xml_file) > 0, "Tests not found.  It's possible your testid, {} is wrong.".format(testid))
    envxml=(EnvBuild(".",infile=xml_file[0]))
    compiler=envxml.get_value("COMPILER")
    mpilib=envxml.get_value("MPILIB")

    #
    # Retrieve machine name
    #
    xml_file=glob.glob("*"+testid+"/env_case.xml")
    envxml=(EnvCase(".",infile=xml_file[0]))
    machine=envxml.get_value("MACH")

    #
    # Retrieve baseline tag to compare to
    #
    xml_file=glob.glob("*"+testid+"/env_test.xml")
    envxml=(EnvTest(".",infile=xml_file[0]))
    baseline = envxml.get_value("BASELINE_NAME_CMP")

    #
    # Create XML header
    #

    testxml=TestReporter()
    testxml.setup_header(tagname,machine,compiler,mpilib,testroot,testtype,baseline)

    #
    # Create lists on tests based on the testid in the testroot directory.
    #
    test_names=glob.glob("*"+testid)
    #
    # Loop over all tests and parse the test results
    #
    test_status={}
    for test_name in test_names:
        if not os.path.isfile(test_name+"/TestStatus"):
            continue
        test_status['COMMENT']=""
        test_status['BASELINE']='----'
        test_status['MEMCOMP']='----'
        test_status['MEMLEAK']='----'
        test_status['NLCOMP']='----'
        test_status['STATUS']='----'
        test_status['TPUTCOMP']='----'
        #
        # Check to see if TestStatus is present, if not then continue
        # I might want to set the status to fail
        #
        try:
            lines = [line.rstrip('\n') for line in open(test_name+"/TestStatus")]
        except (IOError, OSError):
            test_status['STATUS']="FAIL"
            test_status['COMMENT']="TestStatus missing. "
            continue
        #
        # Loop over each line of TestStatus, and check for different types of failures.
        #
        for line in lines:
            if "NLCOMP" in line:
                test_status['NLCOMP']=line[0:4]
            if "MEMLEAK" in line:
                test_status['MEMLEAK']=line[0:4]
            if "MEMCOMP" in line:
                test_status['MEMCOMP']=line[0:4]
            if "BASELINE" in line:
                test_status['BASELINE']=line[0:4]
            if "TPUTCOMP" in line:
                test_status['TPUTCOMP']=line[0:4]
                if "FAIL PFS" in line:
                    test_status['STATUS']="FAIL"
            if "INIT" in line:
                test_status['INIT']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['STATUS']="SFAIL"
                    test_status['COMMENT']+="INIT fail! "
                    break
            if "CREATE_NEWCASE" in line:
                test_status['CREATE_NEWCASE']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['STATUS']="SFAIL"
                    test_status['COMMENT']+="CREATE_NEWCASE fail! "
                    break
            if "XML" in line:
                test_status['XML']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['STATUS']="SFAIL"
                    test_status['COMMENT']+="XML fail! "
                    break
            if "SETUP" in line:
                test_status['SETUP']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['STATUS']="SFAIL"
                    test_status['COMMENT']+="SETUP fail! "
                    break
            if "SHAREDLIB_BUILD" in line:
                test_status['SHAREDLIB_BUILD']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['STATUS']="CFAIL"
                    test_status['COMMENT']+="SHAREDLIB_BUILD fail! "
                    break
            if "MODEL_BUILD" in line:
                test_status['MODEL_BUILD']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['STATUS']="CFAIL"
                    test_status['COMMENT']+="MODEL_BUILD fail! "
                    break
            if "SUBMIT" in line:
                test_status['STATUS']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['COMMENT']+="SUBMIT fail! "
                    break
            if "RUN" in line:
                test_status['STATUS']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['COMMENT']+="RUN fail! "
                    break
            if "COMPARE_base_rest" in line:
                test_status['STATUS']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['COMMENT']+="Restart fail! "
                    break
            if "COMPARE_base_hybrid" in line:
                test_status['STATUS']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['COMMENT']+="Hybrid fail! "
                    break
            if "COMPARE_base_multiinst" in line:
                test_status['STATUS']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['COMMENT']+="Multi instance fail! "
                    break
            if "COMPARE_base_test" in line:
                test_status['STATUS']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['COMMENT']+="Base test fail! "
                    break
            if "COMPARE_base_single_thread" in line:
                test_status['STATUS']=line[0:4]
                if line[0:4] in ("FAIL","PEND"):
                    test_status['COMMENT']+="Thread test fail! "
                    break

            #
            #  Do not include time comments.  Just a preference to have cleaner comments in the test database
            #
            try:
                if 'time=' not in line and 'GENERATE' not in line:
                    if 'BASELINE' not in line:
                        test_status['COMMENT']+=line.split(' ',3)[3]+' '
                    else:
                        test_status['COMMENT']+=line.split(' ',4)[4]+' '
            except Exception: # Probably want to be more specific here
                pass

        #
        # Fill in the xml with the test results
        #
        testxml.add_result(test_name,test_status)

    return testxml

##############################################################################
def _main_func():
###############################################################################

    testroot, testid, tagname, testtype, dryrun, dumpxml = parse_command_line(sys.argv)

    testxml = get_testreporter_xml(testroot, testid, tagname, testtype)

    #
    # Dump xml to a file.
    #
    if dumpxml:
        GenericXML.write(testxml,outfile="TestRecord.xml")


    #
    # Prompt for username and password, then post the XML string to the test database website
    #
    if not dryrun:
        testxml.push2testdb()

###############################################################################

if __name__ == "__main__":
    _main_func()
