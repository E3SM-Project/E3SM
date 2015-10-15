#!/usr/bin/env python

import os, shutil, fnmatch
import argparse
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-o", "--core", dest="core", help="Core that conatins configurations to clean", metavar="CORE", required=True)
parser.add_argument("-c", "--configuration", dest="configuration", help="Configuration to clean", metavar="CONFIG", required=True)
parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution of configuration to clean", metavar="RES", required=True)

args = parser.parse_args()


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
