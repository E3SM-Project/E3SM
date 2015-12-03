#!/usr/bin/env python
# Phillip J. Wolfram, 11/17/2015

import os
import subprocess
import fnmatch
import re
import argparse
import glob
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-c", "--clean", dest="clean", help="If set, script will clean all configurations before running.", action='store_true')
parser.add_argument("-f", "--config_file", dest="config_file", help="Configuration file for running tests.", metavar="FILE", required=True)
parser.add_argument("-o", "--core", dest="core", help="Core that contains configurations to be run.", metavar="CORE", required=True)
parser.add_argument("--base_path", dest="base_path", help="If set, script will create case directories in base_path rather than the current directory.", metavar="PATH")

args = parser.parse_args()

# get total number of test cases
testcases = subprocess.check_output(['./list_testcases.py'])
totalnum = int(re.findall('[0-9]?[0-9]:',testcases)[-1][:-1])

if args.clean:
    print '================================================================================'
    print 'Cleaning all test cases:',
    for i in xrange(totalnum):
            subprocess.check_output(['./clean_testcase.py', '-n', '%d'%(int(i+1))])
    print ' done.'
    print '================================================================================'

print '================================================================================'
print 'Setting up the test cases:'
print '================================================================================'
paths = []
for i in xrange(totalnum):
    setuplist = ['./setup_testcase.py', '-n', '%d'%(int(i+1)), '-f', args.config_file]
    if args.base_path:
        setuplist.append('--base_path')
        setuplist.append(args.base_path)
    setup = subprocess.check_output(setuplist)
    print setup
    paths.append(re.findall('in (.*)\n',setup)[0])
print '================================================================================'
print 'Finished setting up the test cases.'
print '================================================================================'

# running all the test cases
print 'Running all the %s test cases:'%(args.core)
basepath = os.getcwd()
for i in xrange(totalnum):
    # get the core
    if args.base_path:
        core =re.findall('/.*?/',paths[i].replace(args.base_path,''))[0][1:-1]
    else:
        core =re.findall('/.*?/',paths[i].replace(basepath,''))[0][1:-1]
    # make sure this in the right core to test
    if core == args.core:
        print '================================================================================'
        print 'Running test case in %s:'%(paths[i])
        print '================================================================================'
        # run the test
        os.chdir(paths[i])
        runscript = None
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, '*.py'):
                runscript = file

        if runscript is not None:
            runoutput = subprocess.check_output(['python', runscript])
            print runoutput
        print '================================================================================'
        print 'Finished running test case in %s.'%(paths[i])
        print '================================================================================'
    else:
        print '================================================================================'
        print 'Skipping test case in %s:'%(paths[i])
        print '================================================================================'
    os.chdir(basepath)

