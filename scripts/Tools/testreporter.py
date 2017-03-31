#!/usr/bin/env python

"""

Simple script to populate CESM test database with test results.  This can be 
run with out any knowledge of CIME.

"""

import os
import glob
import xml.etree.ElementTree as ET
from xml.dom import minidom
import HTMLParser 
import urllib
import argparse



# Parse command line options

parser = argparse.ArgumentParser(description='Arguements for testreporter')
parser.add_argument("--tagname", 
                        help="Name of the tag being tested.")
parser.add_argument("--testid", 
                        help="Test id, ie c2_0_a6g,c2_0_b6g.")
parser.add_argument("--testroot", 
                        help="Root directory for tests to populate the database.")
parser.add_argument("--testtype", 
                        help="Type of test, prealpha or prebeta.")
parser.add_argument("--dryrun",action="store_true", 
                        help="Do a dry run, database will not be populated.")
parser.add_argument("--dumpxml",action="store_true", 
                        help="Dump XML test results to sceen.")
args = parser.parse_args()

# Fill in values needed from xml files (COMPILER, MPILIB, MACH, BASELINE_NAME_CMP.

os.chdir(args.testroot)

xml_file=glob.glob("*"+args.testid+"/env_build.xml")
root = ET.parse(xml_file[0]).getroot()
for child in root.iter('group'):
   for entry in child.iter('entry'):
     id = entry.get('id')
     value = entry.get('value')
     if id == "COMPILER":
        compiler=value
     if id == "MPILIB":
        mpilib = value

xml_file=glob.glob("*"+args.testid+"/env_case.xml")
root = ET.parse(xml_file[0]).getroot()
for child in root.iter('group'):
   for entry in child.iter('entry'):
     id = entry.get('id')
     value = entry.get('value')
     if id == "MACH":
        machine=value

xml_file=glob.glob("*"+args.testid+"/env_test.xml")
root = ET.parse(xml_file[0]).getroot()
for child in root.iter('group'):
   for entry in child.iter('entry'):
     id = entry.get('id')
     value = entry.get('value')
     if id == "BASELINE_NAME_CMP":
        baseline=value


#
# Create XML
#

testrecord    = ET.Element("testrecord")
tag_name      = ET.SubElement(testrecord,'tag_name').text=args.tagname
mach          = ET.SubElement(testrecord,'mach').text=machine
compiler      = ET.SubElement(testrecord,'compiler',attrib={"version":""}).text=compiler
mpilib        = ET.SubElement(testrecord,'mpilib',attrib={"version":""}).text=mpilib
testroot      = ET.SubElement(testrecord,'testroot').text=args.testroot
testtype      = ET.SubElement(testrecord,'testtype').text=args.testtype
baselinetag  = ET.SubElement(testrecord,'baselinetag').text= baseline
#
# Create lists on tests based on the testid in the testroot directory.
#
test_names=glob.glob("*"+args.testid)
#
# Loop over all tests and parse the test results
#
test_status={}
for test_name in test_names:  
    if test_name == "cs.status."+args.testid:
        continue
    test_status[test_name,'COMMENT']=""
    test_status[test_name,'BASELINE']='----'
    test_status[test_name,'MEMCOMP']='----'
    test_status[test_name,'MEMLEAK']='----'
    test_status[test_name,'NLCOMP']='----'
    test_status[test_name,'STATUS']='----'
    test_status[test_name,'TPUTCOMP']='----'
    #
    # Check to see if TestStatus is present, if not then continue
    # I might want to set the status to fail
    #
    try:
        lines = [line.rstrip('\n') for line in open(test_name+"/TestStatus")]
    except:
        test_status[test_name,'STATUS']="FAIL"
        test_status[test_name,'COMMENT']="TestStatus missing. "
        continue
    #
    # Loop over each line of TestStatus, and check for different types of failures.
    #
    for line in lines:
        if "NLCOMP" in line:
           test_status[test_name,'NLCOMP']=line[0:4]
        if "MEMLEAK" in line:
           test_status[test_name,'MEMLEAK']=line[0:4]
        if "MEMCOMP" in line:
           test_status[test_name,'MEMCOMP']=line[0:4]
        if "BASELINE" in line:
           test_status[test_name,'BASELINE']=line[0:4]
        if "TPUTCOMP" in line:
           test_status[test_name,'TPUTCOMP']=line[0:4]
        if "INIT" in line:
           test_status[test_name,'INIT']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'STATUS']="SFAIL"
               test_status[test_name,'COMMENT']+="INIT fail! "
               break
        if "CREATE_NEWCASE" in line:
           test_status[test_name,'CREATE_NEWCASE']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'STATUS']="SFAIL"
               test_status[test_name,'COMMENT']+="CREATE_NEWCASE fail! "
               break
        if "XML" in line:
           test_status[test_name,'XML']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'STATUS']="SFAIL"
               test_status[test_name,'COMMENT']+="XML fail! "
               break
        if "SETUP" in line:
           test_status[test_name,'SETUP']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'STATUS']="SFAIL"
               test_status[test_name,'COMMENT']+="SETUP fail! "
               break
        if "SHAREDLIB_BUILD" in line:
           test_status[test_name,'SHAREDLIB_BUILD']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'STATUS']="CFAIL"
               test_status[test_name,'COMMENT']+="SHAREDLIB_BUILD fail! "
               break
        if "MODEL_BUILD" in line:
           test_status[test_name,'MODEL_BUILD']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'STATUS']="CFAIL"
               test_status[test_name,'COMMENT']+="MODEL_BUILD fail! "
               break
        if "SUBMIT" in line:
           test_status[test_name,'STATUS']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'COMMENT']+="SUBMIT fail! "
               break
        if "RUN" in line:
           test_status[test_name,'STATUS']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'COMMENT']+="RUN fail! "
               break
        if "COMPARE_base_rest" in line:
           test_status[test_name,'STATUS']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'COMMENT']+="Restart fail! "
               break
        if "COMPARE_base_hybrid" in line:
           test_status[test_name,'STATUS']=line[0:4]
           if line[0:4] == "FAIL":
               test_status[test_name,'COMMENT']+="Hybrid fail! "
               break

        #
        #  Do not include time comments.  Just a preference to have cleaner comments in the test database
        #
        try:
           if 'time=' not in line:
               test_status[test_name,'COMMENT']+=line.split(' ',3)[3]+' '
        except:
           pass

    #
    # File in the xml with the test results
    #
    tests      = ET.Element('tests',attrib={"testname":test_name})
    testrecord.append(tests)
    category=ET.SubElement(tests,'category',attrib={"name":"casestatus"})
    category=ET.SubElement(tests,'category',attrib={"name":"comment"}).text= test_status[test_name,'COMMENT']
    category=ET.SubElement(tests,'category',attrib={"name":"compare"}).text= test_status[test_name,'BASELINE']
    category=ET.SubElement(tests,'category',attrib={"name":"memcomp"}).text= test_status[test_name,'MEMCOMP']
    category=ET.SubElement(tests,'category',attrib={"name":"memleak"}).text=test_status[test_name,'MEMLEAK']
    category=ET.SubElement(tests,'category',attrib={"name":"nlcomp"}).text= test_status[test_name,'NLCOMP']
    category=ET.SubElement(tests,'category',attrib={"name":"status"}).text= test_status[test_name,'STATUS']
    category=ET.SubElement(tests,'category',attrib={"name":"tputcomp"}).text= test_status[test_name,'TPUTCOMP']


#
# Convert XML to a string
#
xmlstr = ET.tostring(testrecord,method="xml",encoding="UTF-8")

#
# Make the XML string human readable and print it out
#
xml=minidom.parseString(xmlstr)
testXML = xml.toprettyxml(encoding="UTF-8")
#
# Dump xml to the screen.
#
if args.dumpxml:
  print testXML

#
# Prompt for username and password, then post the XML string to the test database website
#
if not args.dryrun:
    username=raw_input("Username:")
    os.system("stty -echo")
    password=raw_input("Password:")
    os.system("stty echo")
    params={'username':username,'password':password,'testXML':testXML}
    url="https://csegweb.cgd.ucar.edu/testdb/cgi-bin/processXMLtest.cgi"
    params = urllib.urlencode(params)
    f = urllib.urlopen(url, params)
    #
    # Print any messages from the post command
    #
    print f.read()
    print f.code


