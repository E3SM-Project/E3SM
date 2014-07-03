#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections

# parses both nonlinear and linear solver iteration data from output file, 
# creates an ascii file for them to be used by a subsequent plotting routine
# returns iteration counts and filenames where the same data is stored in /data

def jobprocess(file,job_name):    # Reading a job Output File for production run

#Initialize lists
	proclist = []
	procttl = []
	nonlist = []                                     
	list = []
	avg = []
	avg2 = []
	total = 0.000
	average = 0.000
	out_flag = 0

#Open Output File (From job) if there is one inputted
	if file:
		try:
        		logfile = open(file, "r")
		except:
	        	print 'error reading ' + file
	        	sys.exit(1)
			raise

		for line in logfile:

			#calculate total number of processors
       			if ('total procs = ' in line):
       		         		proclist.append(int(line.split()[7]))
       		         		proctotal = sum(proclist)
       	         			procttl.append(proctotal)
	
			#create nonlinear iterations list
       	 		if ('Nonlinear Solver Step' in line):
       	         		current_step = int(line.split()[4])
	
			#create linear iterations list
			if ('The Belos solver' in line):
	       	         	list.append(int(line.split()[21]))
	
			#calculate average number of linear iterations if time step converges
       	 		if ('Converged!' in line):
       	         		nonlist.append(current_step)
       	         		for value in list:
       	               		  	total += value
				average = total / len(list)
				avg.append(average)
				for n in avg:
					avg2.append(str(round(n, 3)))
				list = []
				total = 0.000
				average = 0.000
				avg = []
	
			#if there contains a time step that fails to converge
			elif ('Failed!' in line):
#				nonlist.append(str(current_step) + "***")
       	         		nonlist.append(current_step)
				for value in list:
					total += value
				average = total / len(list)
				avg.append(average)
				for n in avg:
					avg2.append(str(round(n, 3)))
				list = []
				total = 0.000
				average = 0.000
				avg = []
				out_flag = 1
	
			#no information, indicating something else is wrong and needs investigation
			#else:
			#	error_flag = 1
	
	# write nonlinear solver iteration data

        nd_name = '/newton_' + job_name + '.asc'
        ld_name = '/fgmres_' + job_name + '.asc'

        try:
                iter_n = open("data/" + nd_name, "w")
        except:
                print "error reading newton solver iteration count file"
                sys.exit(1)
                raise
        for line in nonlist:
                snonlist = str(line)
                iter_n.write(snonlist +'\n')

        iter_n.close()

	# write linear solver iteration data
        try:
                iter_l = open("data/" + ld_name, "w")
        except:
                print "error reading fgmres solver iteration count file"
                sys.exit(1)
                raise
        for line in avg2:
                savg2 = str(line)
                iter_l.write(savg2 +'\n')

        iter_l.close()

	return procttl, nonlist, avg2, out_flag, nd_name, ld_name


# end of job diagnostic processing

