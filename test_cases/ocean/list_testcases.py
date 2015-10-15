#!/usr/bin/env python

import os, fnmatch
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-l", "--list", dest="list", help="If set, script will list available configurations which can be passed into setup_testcases.py.", action="store_true")

args = parser.parse_args()

print "Available test cases are:"
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
							print "  -o %s -c %s -r %s"%(core_dir, config_dir, res_dir)
