#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_outprocess
import VV_utilities
import VV_checks

def circdetails(solver_file,job_path):  # using data, fill the web page with info

        failedt_list = []
	solver_file.write('<HTML>\n')
	solver_file.write('<H3>Circular Shelf Iteration Count Details:</H3>')
	solver_file.write('<H4>Eventually published in plot form</H4>')
	solver_file.write('<BR> \n')

# JFNK 2 proc

# Failure checking
        failedt1 = VV_checks.failcheck(job_path, '/circular-shelf/data/circular-shelf.out')
        failedt_list.append(failedt1)

        solver_file.write('<H4>New Run: circular-shelf.out</H4>')
        procttl_circd, nonlist_circd, avg2_circd, out_flag_circd, ndcirc_name, ldcirc_name = VV_outprocess.jobprocess(job_path + '/circular-shelf/data/circular-shelf.out','circ')

        solver_file.write("Number of Processors = " + str(procttl_circd[-1]) + "<BR>\n")
        solver_file.write("Number of Nonlinear Iterations = ")
        VV_utilities.format(solver_file, nonlist_circd)
        solver_file.write('<BR>\n')
        if out_flag_circd == 1:
                solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        solver_file.write("Average Number of Linear Iterations per Time-Step = ")
        VV_utilities.format(solver_file, avg2_circd)
        solver_file.write('<BR> \n')

        solver_file.write('<H4>Benchmark Run: circular-shelf.out</H4>')
        procttl_circb, nonlist_circb, avg2_circb, out_flag_circb, ndcircb_name, ldcircb_name = VV_outprocess.jobprocess(job_path + '/bench/circular-shelf/data/circular-shelf.out','circb')

        solver_file.write("Number of Processors = " + str(procttl_circb[-1]) + "<BR>\n")
        solver_file.write("Number of Nonlinear Iterations = ")
        VV_utilities.format(solver_file, nonlist_circb)
        solver_file.write('<BR>\n')
        if out_flag_circb == 1:
                solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        solver_file.write("Average Number of Linear Iterations per Time-Step = ")
        VV_utilities.format(solver_file, avg2_circb)
        solver_file.write('<BR> \n')

        if 1 in failedt_list:
            failedt = 1
        else:
            failedt = 0

        return failedt
# TODO have jobprocess grab picard solver info as well

def confdetails(solver_file,job_path):  # using data, fill the web page with info

        failedt_list = []

	solver_file.write('<HTML>\n')
	solver_file.write('<H3>Confined Shelf Iteration Count Details:</H3>')
	solver_file.write('<H4>Eventually published in plot form</H4>')
	solver_file.write('<BR> \n')

# JFNK 2 proc
        
# Failure checking
        failedt1 = VV_checks.failcheck(job_path, '/confined-shelf/data/confined-shelf.out')
        failedt_list.append(failedt1)

        solver_file.write('<H4>New Run: confined-shelf.out</H4>')
        procttl_confd, nonlist_confd, avg2_confd, out_flag_confd, ndconf_name, ldconf_name = VV_outprocess.jobprocess(job_path + '/confined-shelf/data/confined-shelf.out', 'conf')

        solver_file.write("Number of Processors = " + str(procttl_confd[-1]) + "<BR>\n")
        solver_file.write("Number of Nonlinear Iterations = ")
        VV_utilities.format(solver_file, nonlist_confd)
        solver_file.write('<BR>\n')
        if out_flag_confd == 1:
                solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        solver_file.write("Average Number of Linear Iterations per Time-Step = ")
        VV_utilities.format(solver_file, avg2_confd)
        solver_file.write('<BR> \n')

        solver_file.write('<H4>Benchmark Run: confined-shelf.out</H4>')
        procttl_confb, nonlist_confb, avg2_confb, out_flag_confb, ndconfb_name, ldconfb_name = VV_outprocess.jobprocess(job_path + '/bench/confined-shelf/data/confined-shelf.out','confb')

        solver_file.write("Number of Processors = " + str(procttl_confb[-1]) + "<BR>\n")
        solver_file.write("Number of Nonlinear Iterations = ")
        VV_utilities.format(solver_file, nonlist_confb)
        solver_file.write('<BR>\n')
        if out_flag_confb == 1:
                solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        solver_file.write("Average Number of Linear Iterations per Time-Step = ")
        VV_utilities.format(solver_file, avg2_confb)
        solver_file.write('<BR> \n')

	solver_file.write('</HTML>\n')
	solver_file.close()

        if 1 in failedt_list:
            failedt = 1
        else:
            failedt = 0
        
        return failedt

def circplot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

        tmpath = job_path + '/circular-shelf/data/circular-shelf.gnu.JFNK.nc'
        if VV_utilities.emptycheck(tmpath) == 0:
                return

	plot_file.write('<HTML>\n')
	plot_file.write('<H3>Circular Shelf Plot Details:</H3>')
        circ_plotfile=''+ ncl_path + '/circshelf.ncl'
	stock='STOCK = addfile(\"'+ job_path + '/bench/circular-shelf/data/circular-shelf.gnu.JFNK.nc\", \"r\")'
	VAR1  ='VAR1 = addfile(\"' + job_path + '/circular-shelf/data/circular-shelf.gnu.JFNK.nc\", \"r\")'
	png  = 'PNG = "' + ncl_path + '/circular-shelf"'
        plot_circ = "ncl '" + VAR1 + "' '" + stock + "' '" + png + "' " + circ_plotfile 

#TODO create an iteration plot and have that also in the html file 
        try:
                output = subprocess.call(plot_circ, shell=True)
        except:
                print "error creating ncl circular shelf plot"
                raise

# transferring circ pic to www file

        if (ncl_path + '/circular-shelf.png'):
        	circpic = "mv -f " + ncl_path + "/circular-shelf.png" + " " + html_path + "/"
        	try:
                	output = subprocess.call(circpic, shell=True)
        	except:
                	print "error moving circular shelf png file to www directory"
                        sys.exit(1)
                	raise

        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>Circular Shelf </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark for a range of processor counts for a range of variables</H4>\n')
        plot_file.write('<OBJECT data="circular-shelf.png" type="image/png" width="1100" height="800" hspace=10 align=left alt="Circular Shelf Plots PNG">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

	plot_file.write('</HTML>\n')
	plot_file.close()

def confplot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

        tmpath = job_path + '/confined-shelf/data/confined-shelf.gnu.JFNK.nc'
        if VV_utilities.emptycheck(tmpath) == 0:
                return

        plot_file.write('<HTML>\n')
	plot_file.write('<H3>Confined Shelf Plot Details:</H3>')
	conf_plotfile=''+ ncl_path + '/confshelf.ncl'
	stock='STOCK = addfile(\"'+ job_path + '/bench/confined-shelf/data/confined-shelf.gnu.JFNK.nc\", \"r\")'
	VAR1  ='VAR1 = addfile(\"' + job_path + '/confined-shelf/data/confined-shelf.gnu.JFNK.nc\", \"r\")'
	png  = 'PNG = "' + ncl_path + '/confined-shelf"'
        plot_conf = "ncl '" + VAR1 + "' '" + stock + "' '" + png + "' " + conf_plotfile 

        try:
                output = subprocess.call(plot_conf, shell=True)
        except:
                print "error creating ncl confined shelf plot"
                raise

# transferring conf pic to www file

        if (ncl_path + '/confined-shelf.png'):
        	confpic = "mv -f " + ncl_path + "/confined-shelf.png" + " " + html_path + "/"
        	try:
                	output = subprocess.call(confpic, shell=True)
        	except:
                	print "error moving confined shelf png file to www directory"
                        sys.exit(1)
                	raise

        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>Confined Shelf </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark for a range of processor counts for a range of variables</H4>\n')
        plot_file.write('<OBJECT data="confined-shelf.png" type="image/png" width="1100" height="800" hspace=10 align=left alt="Confined Shelf Plots PNG">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

	plot_file.write('</HTML>\n')
	plot_file.close()
