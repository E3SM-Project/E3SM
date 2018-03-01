#!/usr/bin/env python
"""
This script is used to manage regression suites. A regression suite is a set of
test cases that ensure one or more features in a model meet certain criteria.

Using this script one can setup or clean a regression suite.

When setting up a regression suite, this script will generate a script to run
all tests in the suite, and additionally setup each individual test case.

When cleaning a regression suite, this script will remove any generated files
for each individual test case, and the run script that runs all test cases.
"""

import sys, os, glob, shutil, numpy, math
import fnmatch
import argparse
import xml.etree.ElementTree as ET
import subprocess
import ConfigParser

def process_test_setup(test_tag, config_file, work_dir, model_runtime, suite_script, baseline_dir):#{{{
	dev_null = open('/dev/null', 'a')

	# Process test attributes
	try:
		test_name = test_tag.attrib['name']
	except:
		print "ERROR: <test> tag is missing 'name' attribute."
		print "Exiting..."
		sys.exit(1)
	
	try:
		test_core = test_tag.attrib['core']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'core' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)

	try:
		test_configuration = test_tag.attrib['configuration']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'configuration' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)

	try:
		test_resolution = test_tag.attrib['resolution']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'resolution' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)

	try:
		test_test = test_tag.attrib['test']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'test' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)
	
	# Determine the file name for the test case output
	case_output_name = test_name.replace(' ', '_')
	
	# Setup test case

	if baseline_dir == 'NONE':
		subprocess.check_call(['./setup_testcase.py', '-q', '-f', config_file, '--work_dir', work_dir,
								'-o', test_core, '-c', test_configuration, '-r', test_resolution, '-t', test_test,
								'-m', model_runtime], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
	else:
		subprocess.check_call(['./setup_testcase.py', '-q', '-f', config_file, '--work_dir', work_dir,
								'-o', test_core, '-c', test_configuration, '-r', test_resolution, '-t', test_test,
								'-m', model_runtime, '-b', baseline_dir], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
	
	print "   -- Setup case '%s': -o %s -c %s -r %s -t %s"%(test_name, test_core, test_configuration, test_resolution, test_test)

	# Write step into suite script to cd into the base of the regression suite
	suite_script.write("os.chdir(base_path)\n")

	# Write the step to define the output file
	suite_script.write("case_output = open('case_outputs/%s', 'w')\n"%(case_output_name))

	# Write step to cd into test case directory
	suite_script.write("os.chdir('%s/%s/%s/%s')\n"%(test_core, test_configuration, test_resolution, test_test))

	for script in test_tag:
		# Process test case script
		if script.tag == 'script':
			try:
				script_name = script.attrib['name']
			except:
				print "ERROR: <script> tag is missing 'name' attribute."
				print 'Exiting...'
				sys.exit(1)

			command = "subprocess.check_call(['%s/%s/%s/%s/%s/%s']"%(work_dir, test_core, test_configuration, test_resolution, test_test, script_name)
			command = '%s, stdout=case_output, stderr=case_output'%(command)
			command = '%s, env=os.environ.copy())'%(command)

			# Write test case run step
			suite_script.write("print ' ** Running case %s'\n"%(test_name))
			suite_script.write('try:\n')
			suite_script.write('\t%s\n'%(command))
			suite_script.write("\tprint '      PASS'\n")
			suite_script.write('except:\n')
			suite_script.write("\tprint '   ** FAIL (See case_outputs/%s for more information)'\n"%(case_output_name))
			suite_script.write("\ttest_failed = True\n")
	
	# Finish writing test case output
	suite_script.write("case_output.close()\n")
	suite_script.write("\n")
	dev_null.close()
#}}}

def process_test_clean(test_tag, work_dir, suite_script):#{{{
	dev_null = open('/dev/null', 'a')

	# Process test attributes
	try:
		test_name = test_tag.attrib['name']
	except:
		print "ERROR: <test> tag is missing 'name' attribute."
		print "Exiting..."
		sys.exit(1)
	
	try:
		test_core = test_tag.attrib['core']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'core' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)

	try:
		test_configuration = test_tag.attrib['configuration']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'configuration' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)

	try:
		test_resolution = test_tag.attrib['resolution']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'resolution' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)

	try:
		test_test = test_tag.attrib['test']
	except:
		print "ERROR: <test> tag with name '%s' is missing 'test' attribute."%(test_name)
		print "Exiting..."
		sys.exit(1)
	

	# Clean test case
	subprocess.check_call(['./clean_testcase.py', '-q', '--work_dir', work_dir, '-o', test_core, '-c', 
							test_configuration, '-r', test_resolution, '-t', test_test], stdout=dev_null, stderr=dev_null)
	
	print "   -- Cleaned case '%s': -o %s -c %s -r %s -t %s"%(test_name, test_core, test_configuration, test_resolution, test_test)
	
	dev_null.close()

#}}}

def setup_suite(suite_tag, work_dir, model_runtime, config_file, baseline_dir):#{{{
	try:
		suite_name = suite_tag.attrib['name']
	except:
		print "ERROR: <regression_suite> tag is missing 'name' attribute."
		print 'Exiting...'
		sys.exit(1)

	if not os.path.exists('%s'%(work_dir)):
		os.makedirs('%s'%(work_dir))
	
	# Create regression suite run script
	regression_script_name = '%s/%s.py'%(work_dir, suite_name)
	regression_script = open('%s'%(regression_script_name), 'w')

	# Write script header
	regression_script.write('#!/usr/bin/env python\n')
	regression_script.write('\n')
	regression_script.write('### This script was written by manage_regression_suite.py as part of a regression_suite file\n')
	regression_script.write('\n')
	regression_script.write('import sys, os\n')
	regression_script.write('import subprocess\n')
	regression_script.write('\n')
	regression_script.write("os.environ['PYTHONUNBUFFERED'] = '1'\n")
	regression_script.write("test_failed = False\n")
	regression_script.write('\n')
	regression_script.write("if not os.path.exists('case_outputs'):\n")
	regression_script.write("\tos.makedirs('case_outputs')\n")
	regression_script.write('\n')
	regression_script.write("base_path = '%s'\n"%(work_dir))

	for child in suite_tag:
		# Process <test> tags within the test suite
		if child.tag == 'test':
			process_test_setup(child, config_file, work_dir, model_runtime, regression_script, baseline_dir)
	
	regression_script.write("if test_failed:\n")
	regression_script.write("\tsys.exit(1)\n")
	regression_script.write("else:\n")
	regression_script.write("\tsys.exit(0)\n")
	regression_script.close()

	dev_null = open('/dev/null', 'a')
	subprocess.check_call(['chmod', 'a+x', '%s'%(regression_script_name)], stdout=dev_null, stderr=dev_null, env=os.environ.copy())
	dev_null.close()
#}}}

def clean_suite(suite_tag, work_dir):#{{{
	try:
		suite_name = suite_tag.attrib['name']
	except:
		print "ERROR: <regression_suite> tag is missing 'name' attribute."
		print 'Exiting...'
		sys.exit(1)

	# Remove the regression suite script, if it exists
	regression_script = '%s/%s.py'%(work_dir, suite_name)
	if os.path.exists(regression_script):
		os.remove(regression_script)

	for child in suite_tag:
		# Process <test> children within the <regression_suite>
		if child.tag == 'test':
			process_test_clean(child, work_dir, regression_script)
#}}}

def summarize_suite(suite_tag):#{{{

	max_procs = 1
	max_threads = 1
	max_cores = 1

	for child in suite_tag:
		if child.tag == 'test':
			try:
				test_name = child.attrib['name']
			except:
				print "<test> tag is missing a 'name' attribute"
				print "Exiting..."
				sys.exit(1)

			try:
				test_core = child.attrib['core']
			except:
				print "<test> tag named '%s' is missing a 'core' attribute"%(test_name)
				print "Exiting..."
				sys.exit(1)

			try:
				test_configuration = child.attrib['configuration']
			except:
				print "<test> tag named '%s' is missing a 'configuration' attribute"%(test_name)
				print "Exiting..."
				sys.exit(1)

			try:
				test_resolution = child.attrib['resolution']
			except:
				print "<test> tag named '%s' is missing a 'resolution' attribute"%(test_name)
				print "Exiting..."
				sys.exit(1)

			try:
				test_test = child.attrib['test']
			except:
				print "<test> tag named '%s' is missing a 'test' attribute"%(test_name)
				print "Exiting..."
				sys.exit(1)

			test_path = '%s/%s/%s/%s'%(test_core, test_configuration, test_resolution, test_test)
			# Loop over all files in test_path that have the .xml extension.
			for file in os.listdir('%s'%(test_path)):
				if fnmatch.fnmatch(file, '*.xml'):
					# Build full file name
					config_file = '%s/%s'%(test_path, file)

					config_tree = ET.parse(config_file)
					config_root = config_tree.getroot()

					if config_root.tag == 'config':
						for model_run in config_root.iter('model_run'):
							try:
								procs_str = model_run.attrib['procs']
								procs = int(procs_str)
							except:
								procs = 1

							try:
								threads_str = model_run.attrib['threads']
								threads = int(threads_str)
							except:
								threads = 1
							
							cores = threads * procs

							if procs > max_procs:
								max_procs = procs

							if threads > max_threads:
								max_threads = threads

							if cores > max_cores:
								max_cores = cores

					del config_root
					del config_tree

	print "\n"
	print " Summary of test cases:"
	print "      Maximum MPI tasks used: %d"%(max_procs)
	print "      Maximum OpenMP threads used: %d"%(max_threads)
	print "      Maximum Total Cores used: %d"%(max_cores)
#}}}

if __name__ == "__main__":
	# Define and process input arguments
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-t", "--test_suite", dest="test_suite", help="Path to file containing a test suite to setup", metavar="FILE", required=True)
	parser.add_argument("-f", "--config_file", dest="config_file", help="Configuration file for test case setup", metavar="FILE")
	parser.add_argument("-s", "--setup", dest="setup", help="Option to determine if regression suite should be setup or not.", action="store_true")
	parser.add_argument("-c", "--clean", dest="clean", help="Option to determine if regression suite should be cleaned or not.", action="store_true")
	parser.add_argument("-m", "--model_runtime", dest="model_runtime", help="Definition of how to build model run commands on this machine", metavar="FILE")
	parser.add_argument("-b", "--baseline_dir", dest="baseline_dir", help="Location of baseslines that can be compared to", metavar="PATH")
	parser.add_argument("--work_dir", dest="work_dir", help="If set, script will setup the test suite in work_dir rather in this script's location.", metavar="PATH")

	args = parser.parse_args()

	if not args.config_file:
		print "WARNING: Not configuration file specified. Using the default of 'local.config'"
		args.config_file = 'local.config'

	if not os.path.exists(args.config_file):
		parser.error("Configuration file '%s' does not exist. Please create and setup before running again."%(args.config_file))

	if not args.work_dir:
		args.work_dir = os.path.dirname(os.path.realpath(__file__))

	args.work_dir = os.path.abspath(args.work_dir)

	if not args.model_runtime:
		args.model_runtime = '%s/runtime_definitions/mpirun.xml'%(os.path.dirname(os.path.realpath(__file__)))
		print 'WARNING: No model runtime specified. Using the default of %s'%(args.model_runtime)

	if not args.baseline_dir:
		args.baseline_dir = 'NONE'

	if not args.setup and not args.clean:
		print 'WARNING: Neither the setup (-s/--setup) nor the clean (-c/--clean) flags were provided. Script will perform no actions.'

	write_history = False

	# Parse regression_suite file
	suite_tree = ET.parse(args.test_suite)
	suite_root = suite_tree.getroot()

	# If the file was a <regression_suite> file, process it
	if suite_root.tag == 'regression_suite':
		# If cleaning, clean the suite
		if args.clean:
			print "Cleaning Test Cases:"
			clean_suite(suite_root, args.work_dir)
			write_history = True
		# If setting up, set up the suite
		if args.setup:
			print "\n"
			print "Setting Up Test Cases:"
			setup_suite(suite_root, args.work_dir, args.model_runtime, args.config_file, args.baseline_dir)
			summarize_suite(suite_root)
			write_history = True

	# Write the history of this command to the command_history file, for
	# provenance.
	if write_history:
		# Build variables for history output
		old_dir = os.getcwd()
		os.chdir( os.path.dirname( os.path.realpath(__file__) ) )
		git_version = subprocess.check_output(['git', 'describe', '--tags', '--dirty'])
		git_version = git_version.strip('\n')
		os.chdir(old_dir)
		calling_command = ""
		for arg in sys.argv:
			calling_command = "%s%s "%(calling_command, arg)

		history_file_path = '%s/command_history'%(args.work_dir)
		if os.path.exists(history_file_path):
			history_file = open(history_file_path, 'a')
			history_file.write('\n')
		else:
			history_file = open(history_file_path, 'w')

		history_file.write('***********************************************************************\n')
		history_file.write('git_version: %s\n'%(git_version))
		history_file.write('command: %s\n'%(calling_command))
		history_file.write('***********************************************************************\n')
		history_file.close()
