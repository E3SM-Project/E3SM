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
    print('codeMacs.py::Error. pFUnit requires argparse module provided by python version >= 2.7.')
    print('Quitting!'); quit()

# from mods.pre import pre_If
# from mods.pre import pre_Repeat

import mods.pre
from mods.pre.pre_If import *
from mods.pre.pre_Repeat import *

parser = \
  argparse.ArgumentParser( \
                           description='A set of code macros to aid preproccessing and code generation for pfunit research', \
                           usage='%(prog)s -inFile INFILE [options]' \
      )

parser.add_argument('--inFile', help='The input file.')
parser.add_argument('--outFile', help='The output file. (Standard out is the default.)')
args = parser.parse_args()

if __name__ == '__main__' :
    result = pre(inFile=args.inFile, outFile=args.outFile)
    if args.inFile :
        if not args.outFile :
            print result


