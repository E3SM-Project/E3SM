#!/usr/bin/env python
"""
This script is used to clean one or more test cases that have already been
setup.

It will remove directories / driver scripts that were generated as part of
setting up a test case.
"""

import sys, os, shutil, fnmatch, re
import argparse
import subprocess
import xml.etree.ElementTree as ET

if __name__ == "__main__":
	# Define and process input arguments
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-o", "--core", dest="core", help="Core that conatins configurations to clean", metavar="CORE")
	parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration to clean", metavar="CONFIG")
	parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution of configuration to clean", metavar="RES")
	parser.add_argument("-t", "--test", dest="test", help="Test name within a resolution to clean", metavar="TEST")
	parser.add_argument("-n", "--case_number", dest="case_num", help="Case number to clean, as listed from list_testcases.py. Can be a comma delimited list of case numbers.", metavar="NUM")
	parser.add_argument("-q", "--quiet", dest="quiet", help="If set, script will not write a command_history file", action="store_true")
	parser.add_argument("-a", "--all", dest="clean_all", help="Is set, the script will clean all test cases in the work_dir.", action="store_true")
	parser.add_argument("--work_dir", dest="work_dir", help="If set, script will clean case directories in work_dir rather than the current directory.", metavar="PATH")

	args = parser.parse_args()

	if not args.case_num and ( not args.core and not args.configuration and not args.resolution and not args.test) and not args.clean_all:
		print 'Must be run with either the --case_number argument, the --all argument, or all of the core, configuration, resolution, and test arguments.'
		parser.error(' Invalid configuration. Exiting...')

	if args.case_num and args.core and args.configuration and args.resoltuion and args.test and args.clean_all:
		print 'Can only be configured with either --case_number (-n), --all (-a), or all of --core (-o), --configuration (-c), --resolution (-r), and --test (-t).'
		parser.error(' Invalid configuration. Too many options used. Exiting...')

	if not args.clean_all:
		if args.case_num:
			use_case_list = True
			case_list = args.case_num.split(',')
		else:
			use_case_list = False
			case_list = list()
			case_list.append(0)
	else:
		use_case_list = True
		valid_case = 1
		case_num = 1
		case_list = list()

		regex = re.compile('(\d):')
		core_configuration = subprocess.check_output(['./list_testcases.py'])
		for line in core_configuration.split('\n'):
			if not regex.search(line) == None:
				conf_arr = line.replace(":", " ").split()
				case_num = int(conf_arr[0])
				del conf_arr
				case_list.append(case_num)
		del core_configuration
		del regex

	if not args.work_dir:
		args.work_dir = os.getcwd()

	args.work_dir = os.path.abspath(args.work_dir)

	# Build variables for history output
	old_dir = os.getcwd()
	os.chdir( os.path.dirname( os.path.realpath(__file__) ) )
	git_version = subprocess.check_output(['git', 'describe', '--tags', '--dirty'])
	git_version = git_version.strip('\n')
	os.chdir(old_dir)
	calling_command = ""
	write_history = False
	for arg in sys.argv:
		calling_command = "%s%s "%(calling_command, arg)

	# Iterate over all cases in the case_list.
	# There is only one if the (-o, -c, -r) options were used in place of (-n)
	for case_num in case_list:
		# If we're using a case_list, determine the core, configuration, and
		# resolution for the current test case.
		if use_case_list:
			core_configuration = subprocess.check_output(['./list_testcases.py', '-n', '%d'%(int(case_num))])
			config_options = core_configuration.strip('\n').split(' ')
			args.core = config_options[1]
			args.configuration = config_options[3]
			args.resolution = config_options[5]
			args.test = config_options[7]

		# Setup each xml file in the configuration directory:
		test_path = '%s/%s/%s/%s'%(args.core, args.configuration, args.resolution, args.test)
		work_dir = '%s/%s'%(args.work_dir, test_path)

		# Only write history if we did something...
		write_history = False

		# Loop over all files in test_path that have the .xml extension.
		for file in os.listdir('%s'%(test_path)):
			if fnmatch.fnmatch(file, '*.xml'):
				# Build full file name
				config_file = '%s/%s'%(test_path, file)

				# Parse file
				config_tree = ET.parse(config_file)
				config_root = config_tree.getroot()

				# Process <config> files
				if config_root.tag == 'config':
					case_dir = config_root.attrib['case']

					case_paths = case_dir.split('/')
					# Determine the base directory in the case path, to delete
					case_base = case_paths[0]

					# Delete the top level directory that was created, if it exists.
					if os.path.exists('%s/%s'%(work_dir, case_base)):
						if os.path.isdir('%s/%s'%(work_dir, case_base)):
							shutil.rmtree('%s/%s'%(work_dir, case_base))
							write_history = True
							print ' -- Removed case %s/%s'%(work_dir, case_base)

				# Process <driver_script> files
				elif config_root.tag == 'driver_script':
					script_name = config_root.attrib['name']

					# Delete script if it exists
					if os.path.exists('%s/%s'%(work_dir, script_name)):
						os.remove('%s/%s'%(work_dir, script_name))
						write_history = True
						print ' -- Removed driver script %s/%s'%(work_dir, script_name)

				del config_tree
				del config_root

	# Write the history of this command to the command_history file, for
	# provenance.
	if write_history and not args.quiet:
		history_file_path = '%s/command_history'%(args.work_dir)
		if os.path.exists(history_file_path):
			history_file = open(history_file_path, 'a')
			history_file.write('\n')
		else:
			history_file = open(history_file_path, 'w')

		history_file.write('***********************************************************************\n')
		history_file.write('git_version: %s\n'%(git_version))
		history_file.write('command: %s\n'%(calling_command))
		history_file.write('setup the following cases:\n')
		if use_case_list:
			for case_num in case_list:
				core_configuration = subprocess.check_output(['./list_testcases.py', '-n', '%d'%(int(case_num))])
				config_options = core_configuration.strip('\n').split(' ')
				history_file.write('\n')
				history_file.write('\tcore: %s\n'%(config_options[1]))
				history_file.write('\tconfiguration: %s\n'%(config_options[3]))
				history_file.write('\tresolution: %s\n'%(config_options[5]))
				history_file.write('\ttest: %s\n'%(config_options[7]))
		else:
			history_file.write('core: %s\n'%(args.core))
			history_file.write('configuration: %s\n'%(args.configuration))
			history_file.write('resolution: %s\n'%(args.resolution))
			history_file.write('test: %s\n'%(args.test))

		history_file.write('***********************************************************************\n')
		history_file.close()

