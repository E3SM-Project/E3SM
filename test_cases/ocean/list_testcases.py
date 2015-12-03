#!/usr/bin/env python

import os, fnmatch
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-o", "--core", dest="core", help="Core to search for configurations within", metavar="CORE")
parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration name to search for", metavar="CONFIG")
parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution to search for", metavar="RES")
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

case_num = 1
for core_dir in os.listdir('.'):
	if os.path.isdir(core_dir) and not core_dir == '.git':
		for config_dir in os.listdir(core_dir):
			config_path = '%s/%s'%(core_dir, config_dir)
			if os.path.isdir(config_path):
				for res_dir in os.listdir(config_path):
					res_path = '%s/%s'%(config_path, res_dir)
					if os.path.isdir(res_path):
						print_case = False
						for case_file in os.listdir(res_path):
							if fnmatch.fnmatch(case_file, '*.xml'):
								print_case = True

						if print_case:
							if not quiet:
								if (not args.core) or args.core == core_dir:
									if (not args.configuration) or args.configuration == config_dir:
										if (not args.resolution) or args.resolution == res_dir:
											print "  %d: -o %s -c %s -r %s"%(case_num, core_dir, config_dir, res_dir)
							if quiet and case_num == print_num:
								print "-o %s -c %s -r %s"%(core_dir, config_dir, res_dir)
							case_num = case_num + 1
