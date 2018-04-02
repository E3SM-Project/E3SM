#!/usr/bin/env python

from testing_utils import colour_init, print_colour, final_summary
import argparse
import sys
import os
import xml.etree.ElementTree as ET
import imp

# tests defined
tests = [{"name":"regression"     , "needsBase":True,  "description":"Tests whether development and base MPAS models are bit reproducible."},
         {"name":"restartability" , "needsBase":False, "description":"Tests if restarting the model is bit reproducible."},
         {"name":"parallelism"    , "needsBase":False, "description":"Tests whether different processor numbers is bit reproducible."}]

colour_init()

# command line arguments
parser = argparse.ArgumentParser(description='Test MPAS-Seaice')

parser.add_argument("-d", "--dev",        required=True,  dest="mpasDevelopmentDir",                 help="MPAS development directory to test")
parser.add_argument("-b", "--base",       required=False, dest="mpasBaseDir",                        help="MPAS base directory to compare against")
parser.add_argument("-t", "--testsuite",  required=False, dest="testSuite",                          help="Input test suite xml file")
parser.add_argument("-o", "--domainsdir", required=False, dest="domainsDir",                         help="Domains directory")
parser.add_argument("-a", "--avail",      required=False, dest="avail",         action='store_true', help="Print available tests to stdout")
parser.add_argument("-c", "--check",      required=False, dest="check",         action='store_true', help="Check that the testing system is working")

args = parser.parse_args()

# output available tests
if (args.avail):
    print "Available tests: "
    for test in tests:
        print "  %s: %s" %(test["name"], test["description"])
    sys.exit()

# check dev directory and model exists
if (not os.path.exists(args.mpasDevelopmentDir)):
    print "Requested development MPAS directory does not exist"
    sys.exit()

if (not os.path.exists(args.mpasDevelopmentDir + "/seaice_model")):
    print "Requested development MPAS executable does not exist"
    sys.exit()

mpasDevelopmentDir = os.path.abspath(args.mpasDevelopmentDir)

# check base directory and model exists
if (args.mpasBaseDir != None):

    if (not os.path.exists(args.mpasBaseDir)):
        print "Requested base MPAS directory does not exist"
        sys.exit()

    if (not os.path.exists(args.mpasBaseDir + "/seaice_model")):
        print "Requested base MPAS executable does not exist"
        sys.exit()

    mpasBaseDir = os.path.abspath(args.mpasBaseDir)

# test suite
if (args.testSuite == None):
    scriptDirectory = os.path.dirname(os.path.abspath(__file__))
    testSuite = scriptDirectory + "/testsuites/testsuite.standard.xml"
else:
    if (not os.path.exists(args.testSuite)):
        print "Requested test suite does not exist"
        sys.exit()
    testSuite = args.testSuite

# domains directory
if (args.domainsDir == None):
    domainsDir = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
    if (domainsDir == None):
        print "Environment variable MPAS_SEAICE_DOMAINS_DIR must be set if no domains directory specified"
        sys.exit()
else:
    domainsDir = args.domainsDir
    if (not os.path.exists(domainsDir)):
        print "Requested domains directory does not exist"
        sys.exit()



# perform tests
nTests = 0
nFails = 0

# load the testsuite xml document
tree = ET.parse(testSuite)
testsuite = tree.getroot()

print_colour("Testing MPAS-Seaice", "title")

print "Test suite: ", testSuite
print

# loop over configurations
for configuration in testsuite:

    # loop over domain
    for domain in configuration:

        print "Using configuration " + configuration.get('name') + " and domain " + domain.get('name') + "..."

        # loop over tests
        for test in domain:

            # get test
            foundTest = False
            for testAvail in tests:

                # check test is available
                if (testAvail['name'] == test.get('name')):

                    foundTest = True

                    # gather test options
                    options = {}
                    for option in test:
                        options[option.get('name')] = option.get('value')

                    # run test
                    module = imp.load_source(testAvail["name"], os.path.dirname(os.path.abspath(__file__)) + "/tests/" + testAvail["name"]+".py")
                    test_function = getattr(module, testAvail["name"])
                    if (testAvail["needsBase"]):
                        failed = test_function(mpasDevelopmentDir, mpasBaseDir, domainsDir, domain.get('name'), configuration.get('name'), options, args.check)
                    else:
                        failed = test_function(mpasDevelopmentDir,              domainsDir, domain.get('name'), configuration.get('name'), options, args.check)

                    nTests = nTests + 1
                    nFails = nFails + failed

            # see if test wasnt available
            if (not foundTest):

                print "Requested test %s not available" %(test.get('name'))
                sys.exit()


# print final summary of tests
print_colour("Summary","title")

final_summary(nTests, nFails)
