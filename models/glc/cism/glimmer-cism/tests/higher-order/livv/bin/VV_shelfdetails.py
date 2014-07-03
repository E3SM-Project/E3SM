#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_outprocess

def circdetails(circ_file,job_path):  # using data, fill the web page with info

	circ_file.write('<HTML>\n')
	circ_file.write('<H3>Circular Shelf Iteration Count Details:</H3>')
	circ_file.write('<H4>Eventually published in plot form</H4>')
	circ_file.write('<BR> \n')

# JFNK 2 proc
	circ_file.write('<H4>New Run:</H4>')
        procttl_circd, nonlist_circd,avg2_circd,out_flag_circd,ndcirc_name,ldcirc_name = VV_outprocess.jobprocess(job_path + '/circular-shelf/data/circular-shelf.out','circ')

        circ_file.write("Number of Processors = " + str(procttl_circd[-1]) + "<BR>\n")
        circ_file.write("Number of Nonlinear Iterations = ")
        for item in nonlist_circd:
                circ_file.write(str(item) + ", ")
        circ_file.write('<BR>\n')
        if out_flag_circd == 1:
                circ_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        circ_file.write("Average Number of Linear Iterations per Time-Step = ")
        for item in avg2_circd:
                circ_file.write(str(item) + ", ")
        circ_file.write('<BR> \n')

        circ_file.write('<H4>Benchmark Run:</H4>')
        procttl_circb, nonlist_circb,avg2_circb,out_flag_circb, ndcircb_name,ldcircb_name = VV_outprocess.jobprocess(job_path + '/circular-shelf/bench/circular-shelf.out','circb')

        circ_file.write("Number of Processors = " + str(procttl_circb[-1]) + "<BR>\n")
        circ_file.write("Number of Nonlinear Iterations = ")
        for item in nonlist_circb:
                circ_file.write(str(item) + ", ")
        circ_file.write('<BR>\n')
        if out_flag_circb == 1:
                circ_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        circ_file.write("Average Number of Linear Iterations per Time-Step = ")
        for item in avg2_circb:
                circ_file.write(str(item) + ", ")
        circ_file.write('<BR> \n')

# TODO have jobprocess grab picard solver info as well

def confdetails(conf_file,job_path):  # using data, fill the web page with info

	conf_file.write('<HTML>\n')
	conf_file.write('<H3>Confined Shelf Iteration Count Details:</H3>')
	conf_file.write('<H4>Eventually published in plot form</H4>')
	conf_file.write('<BR> \n')

# JFNK 2 proc
        conf_file.write('<H4>New Run:</H4>')
        procttl_confd, nonlist_confd,avg2_confd,out_flag_confd, ndconf_name, ldconf_name = VV_outprocess.jobprocess(job_path + '/confined-shelf/data/confined-shelf.out', 'conf')

        conf_file.write("Number of Processors = " + str(procttl_confd[-1]) + "<BR>\n")
        conf_file.write("Number of Nonlinear Iterations = ")
        for item in nonlist_confd:
                conf_file.write(str(item) + ", ")
        conf_file.write('<BR>\n')
        if out_flag_confd == 1:
                conf_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        conf_file.write("Average Number of Linear Iterations per Time-Step = ")
        for item in avg2_confd:
                conf_file.write(str(item) + ", ")
        conf_file.write('<BR> \n')

        conf_file.write('<H4>Benchmark Run:</H4>')
        procttl_confb, nonlist_confb,avg2_confb,out_flag_confb, ndconfb_name, ldconfb_name = VV_outprocess.jobprocess(job_path + '/confined-shelf/bench/confined-shelf.out','confb')

        conf_file.write("Number of Processors = " + str(procttl_confb[-1]) + "<BR>\n")
        conf_file.write("Number of Nonlinear Iterations = ")
        for item in nonlist_confb:
                conf_file.write(str(item) + ", ")
        conf_file.write('<BR>\n')
        if out_flag_confb == 1:
                conf_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        conf_file.write("Average Number of Linear Iterations per Time-Step = ")
        for item in avg2_confb:
                conf_file.write(str(item) + ", ")
        conf_file.write('<BR> \n')

	conf_file.write('</HTML>\n')
	conf_file.close()

def circplot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

	plot_file.write('<HTML>\n')
	plot_file.write('<H3>Circular Shelf Plot Details:</H3>')
	circ_plotfile=''+ ncl_path + '/circshelf.ncl'
	stock='STOCK = addfile(\"'+ job_path + '/circular-shelf/bench/circular-shelf.gnu.JFNK.nc\", \"r\")'
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

	plot_file.write('<HTML>\n')
	plot_file.write('<H3>Confined Shelf Plot Details:</H3>')
	conf_plotfile=''+ ncl_path + '/confshelf.ncl'
	stock='STOCK = addfile(\"'+ job_path + '/confined-shelf/bench/confined-shelf.gnu.JFNK.nc\", \"r\")'
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
