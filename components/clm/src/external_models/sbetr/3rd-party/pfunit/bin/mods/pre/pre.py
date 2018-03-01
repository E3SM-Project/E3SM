#!/usr/bin/env python

# python2 - Deprecated in python 2.7+
import imp
try:
    imp.find_module('argparse')
    found = True
except ImportError:
    found = False

# Preferred for python 2.7+, python 3
# import importlib
# argparse_loader = importlib.find_loader('argparse')
# found = argparse_loader is not None

if found:    
    import argparse
else:
    print('pre.py::Error. pFUnit requires argparse module provided by python version >= 2.7.')
    print('Quitting!'); quit()

#####

from pre_If import *
from pre_Repeat import *

parser = argparse.ArgumentParser(description='A preproccessor for pfunit research')

parser.add_argument('--inFile')
parser.add_argument('--outFile')
args = parser.parse_args()

if __name__ == '__main__' :
    result = pre(inFile=args.inFile, outFile=args.outFile)
    if args.inFile :
        if not args.outFile :
            print result


