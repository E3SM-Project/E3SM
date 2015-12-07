#!/usr/bin/env python

import os, fnmatch
import argparse

# Define and process input arguments
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-o", "--core", dest="core", help="Core to search for configurations within", metavar="CORE")
parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration name to search for", metavar="CONFIG")
parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution to search for", metavar="RES")
parser.add_argument("-t", "--test", dest="test", help="Test name to search for", metavar="TEST")
parser.add_argument("-n", "--number", dest="number", help="If set, script will print the flags to use a the N'th configuraiton.")

args = parser.parse_args()

quiet = False

try:
	print_num = 0
	if args.number:
		quiet = True
		print_num = int(args.number)
except:
	args.number = 0
	print_num = 0

if not quiet:
	print "Available test cases are:"

# Start case numbering at 1
case_num = 1

# Iterate over all cores
for core_dir in os.listdir('.'):
	if os.path.isdir(core_dir) and not core_dir == '.git':
		# Iterate over all configurations within a core
		for config_dir in os.listdir(core_dir):
			config_path = '%s/%s'%(core_dir, config_dir)
			if os.path.isdir(config_path):
				# Iterate over all resolutions within a configuration
				for res_dir in os.listdir(config_path):
					res_path = '%s/%s'%(config_path, res_dir)
					if os.path.isdir(res_path):
						for test_dir in os.listdir(res_path):
							test_path = '%s/%s'%(res_path, test_dir)
							if os.path.isdir(test_path):
								print_case = False
								# Iterate over all files within a resolution
								for case_file in os.listdir(test_path):
									if fnmatch.fnmatch(case_file, '*.xml'):
										print_case = True

								# Print the options if a case file was found.
								if print_case:
									if not quiet:
										if (not args.core) or args.core == core_dir:
											if (not args.configuration) or args.configuration == config_dir:
												if (not args.resolution) or args.resolution == res_dir:
													if (not args.test) or args.test == test_dir:
														print "  %d: -o %s -c %s -r %s -t %s"%(case_num, core_dir, config_dir, res_dir, test_dir)
									if quiet and case_num == print_num:
										print "-o %s -c %s -r %s -t %s"%(core_dir, config_dir, res_dir, test_dir)
									case_num = case_num + 1
