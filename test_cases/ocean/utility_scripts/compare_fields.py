#!/usr/bin/env python
import sys, os, glob, shutil
import numpy as np

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-1", "--file1", dest="filename1", help="first input file", metavar="FILE")
parser.add_argument("-2", "--file2", dest="filename2", help="second input file", metavar="FILE")
parser.add_argument("-v", "--var", dest="variable", help="variable to compute error with", metavar="VAR")
parser.add_argument("--l2", dest="l2_norm", help="value of L2 norm for a pass.", metavar="VAL")
parser.add_argument("--l1", dest="l1_norm", help="value of L1 norm for a pass.", metavar="VAL")
parser.add_argument("--linf", dest="linf_norm", help="value of L_Infinity norm for a pass.", metavar="VAL")
parser.add_argument("-q", "--quiet", dest="quiet", help="turns off printing if diff passes test.", action="store_true")

args = parser.parse_args()

if not args.filename1:
	parser.error("Two filenames are required inputs.")

if not args.filename2:
	parser.error("Two filenames are required inputs.")

if not args.variable:
	parser.error("Variable is a required input.")

if (not args.l2_norm) and (not args.l1_norm) and (not args.linf_norm):
	print "WARNING: Script will pass since no norm values have been defined."

f1 = NetCDFFile(args.filename1,'r')
f2 = NetCDFFile(args.filename2,'r')

try:
	time_length = f1.variables['xtime'].shape[0]
except:
	time_length = 1

try:
	field1 = f1.variables[args.variable]
	field2 = f2.variables[args.variable]
except:
	print "ERROR: Field '%s' does not exist in both"%(args.variable)
	print "           file1: %s"%(args.filename1)
	print "       and file2: %s"%(args.filename2)
	print "Exiting with a successful comparision, since no comparision can be done."
	sys.exit(0)

if not field1.shape == field2.shape:
	print "ERROR: Field sizes don't match in different files."

linf_norm = -(sys.float_info.max)

pass_val = True

print "Comparing field '%s'"%(args.variable)

for t in range( 0, time_length):
	pass_time = True
	diff = field1[t][:] - field2[t][:]
	l2_norm = sum(diff * diff)
	l2_norm = l2_norm / np.sum(field1[t][:].shape)
	l2_norm = np.max(l2_norm)
	
	l1_norm = sum(diff)
	l1_norm = l1_norm / np.sum(field1[t][:].shape)
	l1_norm = np.max(l1_norm)

	if np.amax(diff) > linf_norm:
		linf_norm = np.amax(diff)

	diff_str = '%d: '%(t)
	if args.l1_norm:
		if float(args.l1_norm) < abs(l1_norm):
			pass_time = False
	diff_str = '%s l1: %16.14e '%(diff_str, l1_norm)

	if args.l2_norm:
		if float(args.l2_norm) < abs(l2_norm):
			pass_time = False
	diff_str = '%s l2: %16.14e '%(diff_str, l2_norm)

	if args.linf_norm:
		if float(args.linf_norm) < abs(linf_norm):
			pass_time = False
	diff_str = '%s linf: %16.14e '%(diff_str, linf_norm)

	if not args.quiet:
		print diff_str
	elif not pass_time:
		print diff_str
	
	if not pass_time:
		pass_val = False

	del diff

if pass_val:
	sys.exit(0)
else:
	sys.exit(1)
