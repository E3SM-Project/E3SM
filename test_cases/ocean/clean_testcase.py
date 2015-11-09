#!/usr/bin/env python

import sys, os, shutil, fnmatch
import argparse
import subprocess
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-o", "--core", dest="core", help="Core that conatins configurations to clean", metavar="CORE")
parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration to clean", metavar="CONFIG")
parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution of configuration to clean", metavar="RES")
parser.add_argument("-n", "--case_number", dest="case_num", help="Case number to clean, as listed from list_testcases.py.", metavar="NUM")
parser.add_argument("--base_path", dest="base_path", help="If set, script will clean case directories in base_path rather than the current directory.", metavar="PATH")

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

if not args.base_path:
	args.base_path = os.getcwd()

# Build variables for history output
git_version = subprocess.check_output(['git', 'describe', '--tags', '--dirty'])
git_version = git_version.strip('\n')
calling_command = ""
for arg in sys.argv:
	calling_command = "%s%s "%(calling_command, arg)


test_path = '%s/%s/%s'%(args.core, args.configuration, args.resolution)
base_path = '%s/%s'%(args.base_path, test_path)
write_history = False
for file in os.listdir('%s'%(test_path)):
	if fnmatch.fnmatch(file, '*.xml'):
		config_file = '%s/%s'%(test_path, file)

		config_tree = ET.parse(config_file)
		config_root = config_tree.getroot()

		if config_root.tag == 'config':
			case_dir = config_root.attrib['case']

			if os.path.exists('%s/%s'%(base_path, case_dir)):
				if os.path.isdir('%s/%s'%(base_path, case_dir)):
					shutil.rmtree('%s/%s'%(base_path, case_dir))
					write_history = True
					print ' -- Removed case %s/%s'%(base_path, case_dir)

		elif config_root.tag == 'driver_script':
			script_name = config_root.attrib['name']

			if os.path.exists('%s/%s'%(base_path, script_name)):
				os.remove('%s/%s'%(base_path, script_name))
				write_history = True
				print ' -- Removed driver script %s/%s'%(base_path, script_name)

		del config_tree
		del config_root

if write_history:
	history_file_path = '%s/command_history'%(args.base_path)
	if os.path.exists(history_file_path):
		history_file = open(history_file_path, 'a')
		history_file.write('\n')
	else:
		history_file = open(history_file_path, 'w')
	
	history_file.write('git_version: %s\n'%(git_version))
	history_file.write('command: %s\n'%(calling_command))
	history_file.write('core: %s\n'%(args.core))
	history_file.write('configuration: %s\n'%(args.configuration))
	history_file.write('resolution: %s\n'%(args.resolution))
	history_file.close()
