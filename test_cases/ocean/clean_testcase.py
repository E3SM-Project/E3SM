#!/usr/bin/env python

import os, shutil, fnmatch
import argparse
import subprocess
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-o", "--core", dest="core", help="Core that conatins configurations to clean", metavar="CORE")
parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration to clean", metavar="CONFIG")
parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution of configuration to clean", metavar="RES")
parser.add_argument("-n", "--case_number", dest="case_num", help="Case number to clean, as listed from list_testcases.py.", metavar="NUM")

args = parser.parse_args()

if not args.case_num and ( not args.core and not args.configuration and not args.resolution):
	print 'Must be run with either the --case_number argument, or the core, configuration, and resolution arguments.'
	parser.error(' Invalid configuration. Exiting...')

if args.case_num and args.core and args.configuration and args.resoltuion:
	print 'Can only be configured with either --case_number (-n) or --core (-o), --configuration (-c), and --resolution (-r).'
	parser.error(' Invalid configuration. Too many options used. Exiting...')

if args.case_num:
	core_configuration = subprocess.check_output(['./list_testcases.py', '-n', '%d'%(int(args.case_num))])
	config_options = core_configuration.strip('\n').split(' ')
	args.core = config_options[1]
	args.configuration = config_options[3]
	args.resolution = config_options[5]

base_path = '%s/%s/%s'%(args.core, args.configuration, args.resolution)
for file in os.listdir('%s'%(base_path)):
	if fnmatch.fnmatch(file, '*.xml'):
		config_file = '%s/%s'%(base_path, file)

		config_tree = ET.parse(config_file)
		config_root = config_tree.getroot()

		case_dir = config_root.attrib['case']

		del config_tree
		del config_root

		if os.path.exists('%s/%s'%(base_path, case_dir)):
			shutil.rmtree('%s/%s'%(base_path, case_dir))
			print " -- Removed case %s/%s"%(base_path, case_dir)
