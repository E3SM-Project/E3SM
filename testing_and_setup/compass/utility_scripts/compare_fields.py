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

files_exist = True

if not os.path.exists(args.filename1):
	print "ERROR: File %s does not exist. Comparison will FAIL."%(args.filename1)
	files_exist = False

if not os.path.exists(args.filename2):
	print "ERROR: File %s does not exist. Comparison will FAIL."%(args.filename2)
	files_exist = False

if not files_exist:
	sys.exit(1)

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
	print "Exiting with a failed comparision, since no comparision can be done but a comparison was requested."
	sys.exit(1)

if not field1.shape == field2.shape:
	print "ERROR: Field sizes don't match in different files."
	sys.exit(1)

linf_norm = -(sys.float_info.max)

pass_val = True

print "Beginning variable comparisons for all time levels of field '%s'. Note any time levels reported are 0-based."%(args.variable)
if ( args.l2_norm or args.l1_norm or args.linf_norm ):
	print "    Pass thresholds are:"
	if ( args.l1_norm ):
		print "       L1: %16.14e"%(float(args.l1_norm))
	if ( args.l2_norm ):
		print "       L2: %16.14e"%(float(args.l2_norm))
	if ( args.linf_norm ):
		print "       L_Infinity: %16.14e"%(float(args.linf_norm))

field_dims = field1.dimensions

if "Time" in field_dims:
	for t in range( 0, time_length):
		pass_time = True
		if len(field_dims) >= 2:
			diff = np.absolute(field1[t][:] - field2[t][:])
		else:
			diff = np.absolute(field1[t] - field2[t])

		l2_norm = np.sum(diff * diff)
		l2_norm = np.sqrt(l2_norm)

		l1_norm = np.sum(diff)
		if len(field_dims) >= 2:
			l1_norm = l1_norm / np.sum(field1[t][:].shape)
		l1_norm = np.max(l1_norm)

		if np.amax(diff) > linf_norm:
			linf_norm = np.amax(diff)

		diff_str = '%d: '%(t)
		if args.l1_norm:
			if float(args.l1_norm) < l1_norm:
				pass_time = False
		diff_str = '%s l1: %16.14e '%(diff_str, l1_norm)

		if args.l2_norm:
			if float(args.l2_norm) < l2_norm:
				pass_time = False
		diff_str = '%s l2: %16.14e '%(diff_str, l2_norm)

		if args.linf_norm:
			if float(args.linf_norm) < linf_norm:
				pass_time = False
		diff_str = '%s linf: %16.14e '%(diff_str, linf_norm)

		if not args.quiet:
			print diff_str
		elif not pass_time:
			print diff_str

		if not pass_time:
			pass_val = False

		del diff
else:
	if len(field_dims) >= 2:
		diff = np.absolute(field1[:] - field2[:])
	else:
		diff = np.absolute(field1[0] - field2[0])

	l2_norm = np.sum(diff * diff)
	l2_norm = np.sqrt(l2_norm)

	l1_norm = np.sum(diff)
	if len(field_dims) >= 2:
		l1_norm = l1_norm / np.sum(field1[:].shape)
	l1_norm = np.max(l1_norm)

	if np.amax(diff) > linf_norm:
		linf_norm = np.amax(diff)

	diff_str = ''
	if args.l1_norm:
		if float(args.l1_norm) < l1_norm:
			pass_val = False
	diff_str = '%s l1: %16.14e '%(diff_str, l1_norm)

	if args.l2_norm:
		if float(args.l2_norm) < l2_norm:
			pass_val = False
	diff_str = '%s l2: %16.14e '%(diff_str, l2_norm)

	if args.linf_norm:
		if float(args.linf_norm) < linf_norm:
			pass_val = False
	diff_str = '%s linf: %16.14e '%(diff_str, linf_norm)

	if not args.quiet:
		print diff_str
	elif not pass_val:
		print diff_str

	del diff

if pass_val:
	sys.exit(0)
else:
	sys.exit(1)
