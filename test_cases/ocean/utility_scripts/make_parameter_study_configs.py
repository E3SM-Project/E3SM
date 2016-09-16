#!/usr/bin/env python
"""
Writes out a series of config files to be used to perform a parameter
study.  A template file, specified with the -t flag, contains dummy
strings prefixed with '@' to be replaced with parameter values from a list.
The resulting config files are numbered consecutively with a prefix povided
with the -o flag.  Parameter names and values are provided as a list with
the -p flag using syntax as in this example:
-p param1=1,2,3 param2=1e3,1e4,1e5 param3='a','b','c' \\
   param4=.true.,.false.,.true.
The number of parameter values must be the same for all parameters and all
parameters are varied simultaneously.
"""

import argparse

def write_from_template(inFile,outFile,replacements):
    inID = open(inFile)
    outID = open(outFile, 'w')

    for line in inID:
        for src, target in replacements.iteritems():
            line = line.replace(src, target)
        outID.write(line)
    inID.close()
    outID.close()

# Define and process input arguments
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-t", "--template", dest="template", help="A config file in which to add or modify a parameter", metavar="TEMPLATE", required=True)
parser.add_argument("-o", "--out_prefix", dest="out_prefix", help="The prefix for the output config file", metavar="PREFIX", required=True)
parser.add_argument("-p", "--parameters", dest="parameters", help="A list of parameters and comma-separated values", metavar="PARAMETERS", nargs="+", required=True)

args = parser.parse_args()

parameters = {}

first = True
for parameterString in args.parameters:
    (parameter, valueString) = parameterString.split('=',1)
    values = valueString.split(',')
    if first:
        valueCount = len(values)
        first = False
    else:
        assert(len(values) == valueCount)

    parameters[parameter] = values

for valueIndex in range(valueCount):
    outFileName = '%s_%02i.xml'%(args.out_prefix, valueIndex)
    replacements = {}
    for parameter in parameters:
        replacements['@%s'%parameter] = parameters[parameter][valueIndex]
    write_from_template(args.template, outFileName, replacements)
