#!/usr/bin/python

import sys
import os
import re

filepath = "/tmp/work/ab3/higher-order/reg_test/circular-shelf/"
filename= "trilinosOptions.xml"
xml_file = filepath + filename

benchpath = "/tmp/work/ab3/higher-order/reg_test/confined-shelf/"
benchname = "trilinosOptions.xml"
bench_xml_file = benchpath + benchname

file_keywords_l1 = ('Convergence Tolerance', 'Maximum Iterations', 'Preconditioner Type', 'Prec Type', 'Overlap', 'fact: level-of-fill')
file_keywords_l2 = ('Jacobian Operator', 'Forcing Term Method', 'Forcing Term Alpha', 'Forcing Term Gamma', 'Matrix-Free Perturbation', 'Linear Solver Type', 'Solver Type', 'Convergence Tolerance', 'Flexible Gmres')

bench_keywords_l1 = ('Convergence Tolerance', 'Maximum Iterations', 'Preconditioner Type', 'Prec Type', 'Overlap', 'fact: level-of-fill')
bench_keywords_l2 = ('Jacobian Operator', 'Forcing Term Method', 'Forcing Term Alpha', 'Forcing Term Gamma', 'Matrix-Free Perturbation', 'Linear Solver Type', 'Solver Type', 'Convergence Tolerance', 'Flexible Gmres')


#create dictionaries
fileL1parameters = {}
fileL2parameters = {}
benchL1parameters = {}
benchL2parameters = {}
file_max_it = []
bench_max_it = []
flag = False

#open file
try:
	xmlfile = open(xml_file, "r")
except:
	print "error reading xml file"
	sys.exit(1)

#search for keywords_l1 (Preconditioner keywords)
while True:
	line = xmlfile.readline()
	if 'Trilinos nonlinear solver options' in line:
		break;
	for kw1 in file_keywords_l1:
		if ('"' + kw1 + '"' in line):
			value = line.split('value=')
			match = re.search('"(.*)"', value[1])
			fileL1parameters[kw1] = match.group(1)	


print fileL1parameters['Convergence Tolerance']


#print dictionary
print 'XML File:'
print 'Preconditioner: Picard'
#for parameter in fileL1parameters.keys():
#	print parameter
print 'Block GMRES: Convergence Tolerance = ', fileL1parameters['Convergence Tolerance']
print 'Block GMRES: Maximum Iterations = ', fileL1parameters['Maximum Iterations']
print 'Preconditioner Type = ', fileL1parameters['Preconditioner Type']
print 'Prec Type = ', fileL1parameters['Prec Type']
print 'Overlap = ', fileL1parameters['Overlap']
print 'fact: level-of-fill = ', fileL1parameters['fact: level-of-fill']

#search for keywords_l2 (Solver keywords)
for line in xmlfile.readlines():
	for kw2 in file_keywords_l2:
		if ('"' + kw2 + '"' in line):
			value = line.split('value=')
			match = re.search('"(.*)"', value[1])
			fileL2parameters[kw2] = match.group(1)
	if 'Forcing Term Alpha' in line:
		flag = True
	if 'Maximum Iterations' in line:
		value = line.split('value=')
		match = re.search('"(.*)"', value[1])
		file_max_it.append(match.group(1))
	xmlfile.close()

#print dictionary
print ' '
print 'Solver: NK'
#for parameter in L2parameters.keys():
#	print parameter + ' = ' + L2parameters[parameter]
print 'Newton: Jacobian Operator = ', fileL2parameters['Jacobian Operator']
print 'Newton: Forcing Term Method = ', fileL2parameters['Forcing Term Method']
print 'Newton: Maximum Iterations = ', file_max_it[1]
if flag == True:
	print 'Newton: Forcing Term Alpha = ', fileL2parameters['Forcing Term Alpha']	
	print 'Newton: Forcing Term Gamma = ', fileL2parameters['Forcing Term Gamma']
print 'Matrix-Free Perturbation = ', fileL2parameters['Matrix-Free Perturbation']
print 'Linear Solver Type = ', fileL2parameters['Linear Solver Type']
print 'Solver Type = ', fileL2parameters['Solver Type']
print 'GMRES: Convergence Tolerance = ', fileL2parameters['Convergence Tolerance']
print 'GMRES: Maximum Iterations = ', file_max_it[0]
print 'GMRES: Flexible GMRES = ', fileL2parameters['Flexible Gmres']

#reset flag for bench file
flag = False


#open file
try:
        benchxmlfile = open(bench_xml_file, "r")
except:
        print "error reading bench xml file"
        sys.exit(1)

#search for keywords_l1 (Preconditioner keywords)
while True:
        line = benchxmlfile.readline()
        if 'Trilinos nonlinear solver options' in line:
                break;
        for kw1 in bench_keywords_l1:
                if ('"' + kw1 + '"' in line):
                        value = line.split('value=')
                        match = re.search('"(.*)"', value[1])
                        benchL1parameters[kw1] = match.group(1)  

#print bench dictionary
print ' '
print ' '
print ' '
print 'Bench XML File:'
print 'Preconditioner: Picard'
#for parameter in L1parameters.keys():
#       print parameter + ' = ' + L1parameters[parameter]
print 'Block GMRES: Convergence Tolerance = ', benchL1parameters['Convergence Tolerance']
print 'Block GMRES: Maximum Iterations = ', benchL1parameters['Maximum Iterations']
print 'Preconditioner Type = ', benchL1parameters['Preconditioner Type']
print 'Prec Type = ', benchL1parameters['Prec Type']
print 'Overlap = ', benchL1parameters['Overlap']
print 'fact: level-of-fill = ', benchL1parameters['fact: level-of-fill']

#search for keywords_l2 (Solver keywords)
for line in benchxmlfile.readlines():
        for kw2 in bench_keywords_l2:
                if ('"' + kw2 + '"' in line):
                        value = line.split('value=')
                        match = re.search('"(.*)"', value[1])
                        benchL2parameters[kw2] = match.group(1)
        if 'Forcing Term Alpha' in line:
                flag = True
	if 'Maximum Iterations' in line:
                value = line.split('value=')
                match = re.search('"(.*)"', value[1])
                bench_max_it.append(match.group(1))
	benchxmlfile.close()

#print dictionary
print ' '
print 'Solver: NK'
#for parameter in L2parameters.keys():
#       print parameter + ' = ' + L2parameters[parameter]
print 'Newton: Jacobian Operator = ', benchL2parameters['Jacobian Operator']
print 'Newton: Forcing Term Method = ', benchL2parameters['Forcing Term Method']
print 'Newton: Maximum Iterations = ', bench_max_it[1]
if flag == True:
        print 'Newton: Forcing Term Alpha = ', benchL2parameters['Forcing Term Alpha']
        print 'Newton: Forcing Term Gamma = ', benchL2parameters['Forcing Term Gamma']
print 'Matrix-Free Perturbation = ', benchL2parameters['Matrix-Free Perturbation']
print 'Linear Solver Type = ', benchL2parameters['Linear Solver Type']
print 'Solver Type = ', benchL2parameters['Solver Type']
print 'GMRES: Convergence Tolerance = ', benchL2parameters['Convergence Tolerance']
print 'GMRES: Maximum Iterations = ', bench_max_it[0]
print 'GMRES: Flexible GMRES = ', benchL2parameters['Flexible Gmres']

