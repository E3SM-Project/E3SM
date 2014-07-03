#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_outprocess

# routine for ISMIP HOM A 80km
def a80details(solver_file,job_path):  # using data, fill the web page with info

	solver_file.write('<HTML>\n')
	solver_file.write('<H3>ISMIP HOM A 80km Iteration Count Details:</H3>')
	solver_file.write('<H4>Eventually published in plot form</H4>')
	solver_file.write('<BR> \n')
# JFNK gnu 1 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_ih1d, nonlist_ih1d,avg2_ih1d,out_flag_ih1d, ndiha1_name, ldiha1_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-a/80km/data/ishom.a.80km.out.1', 'imhoma1')

	solver_file.write("Number of Processors = " + str(procttl_ih1d[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih1d:
		solver_file.write(str(item) + ", ")
# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR>\n')
	if out_flag_ih1d == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih1d:
		solver_file.write(str(item) + ", ")
# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_ih1b, nonlist_ih1b,avg2_ih1b,out_flag_ih1b,ndiha1b_name, ldiha1b_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-a/80km/bench/ishom.a.80km.out.1', 'imhoma1b')

	solver_file.write("Number of Processors = " + str(procttl_ih1b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih1b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih1b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih1b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

# JFNK gnu 2 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_ih2d, nonlist_ih2d,avg2_ih2d,out_flag_ih2d,ndiha2_name,ldiha2_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-a/80km/data/ishom.a.80km.out.2','imhoma2')

	solver_file.write("Number of Processors = " + str(procttl_ih2d[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih2d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih2d == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih2d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_ih2b, nonlist_ih2b,avg2_ih2b,out_flag_ih2b,ndiha2b_name,ldiha2b_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-a/80km/bench/ishom.a.80km.out.2','imhoma2b')

	solver_file.write("Number of Processors = " + str(procttl_ih2b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih2b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih2b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih2b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

# JFNK gnu 4 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_ih4d, nonlist_ih4d,avg2_ih4d,out_flag_ih4d,ndiha4_name,ldiha4_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-a/80km/data/ishom.a.80km.out.4','imhoma4')

	solver_file.write("Number of Processors = " + str(procttl_ih4d[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih4d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih4d == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih4d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_ih4b, nonlist_ih4b,avg2_ih4b,out_flag_ih4b,ndiha4b_name,ldiha4b_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-a/80km/bench/ishom.a.80km.out.4','imhoma4b')

	solver_file.write("Number of Processors = " + str(procttl_ih4b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih4b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih4b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih4b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('</HTML>\n')
	solver_file.close()

# routine for ISMIP HOM C 80km
def c80details(solver_file,job_path):  # using data, fill the web page with info

	solver_file.write('<HTML>\n')
	solver_file.write('<H3>ISMIP HOM C Iteration Count Details:</H3>')
	solver_file.write('<H4>Eventually published in plot form</H4>')
	solver_file.write('<BR> \n')
# JFNK gnu 1 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_ih1d, nonlist_ih1d,avg2_ih1d,out_flag_ih1d,ndihc1_name,ldihc1_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-c/80km/data/ishom.c.80km.out.1','imhomc1')

	solver_file.write("Number of Processors = " + str(procttl_ih1d[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih1d:
		solver_file.write(str(item) + ", ")
# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR>\n')
	if out_flag_ih1d == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih1d:
		solver_file.write(str(item) + ", ")
# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_ih1b, nonlist_ih1b,avg2_ih1b,out_flag_ih1b,ndihc1b_name,ldihc1b_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-c/80km/bench/ishom.c.80km.out.1','imhomc1b')

	solver_file.write("Number of Processors = " + str(procttl_ih1b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih1b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih1b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih1b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

# JFNK gnu 2 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_ih2d, nonlist_ih2d,avg2_ih2d,out_flag_ih2d,ndihc2_name,ldihc2_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-c/80km/data/ishom.c.80km.out.2','imhom2')

	solver_file.write("Number of Processors = " + str(procttl_ih2d[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih2d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih2d == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih2d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_ih2b, nonlist_ih2b,avg2_ih2b,out_flag_ih2b,ndihc2b_name,ldihc2b_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-c/80km/bench/ishom.c.80km.out.2','imhom2b')

	solver_file.write("Number of Processors = " + str(procttl_ih2b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih2b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih2b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih2b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

# JFNK gnu 4 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_ih4d, nonlist_ih4d,avg2_ih4d,out_flag_ih4d,ndihc4_name,ldihc4_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-c/80km/data/ishom.c.80km.out.4','imhomc4')

	solver_file.write("Number of Processors = " + str(procttl_ih4d[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih4d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih4d == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih4d:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_ih4b, nonlist_ih4b,avg2_ih4b,out_flag_ih4b,ndihc4b_name,ldihc4b_name = VV_outprocess.jobprocess(job_path + '/ismip-hom-c/80km/bench/ishom.c.80km.out.4','imhomc4b')

	solver_file.write("Number of Processors = " + str(procttl_ih4b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_ih4b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_ih4b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_ih4b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('</HTML>\n')
	solver_file.close()

def a80plot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

        plot_file.write('<HTML>\n')
        plot_file.write('<H3>ISMIP HOM A 80km Plot Details:</H3>')
        ishoma_plotfile=''+ ncl_path + '/ismipa80.ncl'
        stock='STOCK = addfile(\"'+ job_path + '/ismip-hom-a/80km/bench/ishom.a.80km.JFNK.out.nc\", \"r\")'
        VAR1  ='VAR1 = addfile(\"' + job_path + '/ismip-hom-a/80km/data/ishom.a.80km.JFNK.out.nc\", \"r\")'
        png  = 'PNG = "' + ncl_path + '/ismipa80"'
        plot_ishoma = "ncl '" + VAR1 + "' '" + stock + "' '" + png + "' " + ishoma_plotfile

#TODO create an iteration plot and have that also in the html file 
        try:
                output = subprocess.call(plot_ishoma, shell=True)
        except:
                print "error creating ncl ismip hom a plots"
                raise

# transferring ismipa pic to www file

        if (ncl_path + '/ismipa80.png'): 
                ishomapic = "mv -f " + ncl_path + "/ismipa80.png" + " " + html_path + "/"
                try:
                        output = subprocess.call(ishomapic, shell=True)
                except:
                        print "error moving ismip hom a 80km png file to www directory"
                        sys.exit(1)
                        raise

        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>ISMIP HOM A 80km </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark run</H4>\n')
        plot_file.write('<BR>\n')
        plot_file.write('<OBJECT data="ismipa80.png" type="image/png" width="1100" height="800" hspace=10 align=left alt="ISMIP HOM A 80km Plots">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

        plot_file.write('</HTML>\n')
        plot_file.close()

def c80plot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

        plot_file.write('<HTML>\n')
        plot_file.write('<H3>ISMIP HOM C 80km Plot Details:</H3>')
        ishomc_plotfile=''+ ncl_path + '/ismipc80.ncl'
        stock='STOCK = addfile(\"'+ job_path + '/ismip-hom-c/80km/bench/ishom.c.80km.JFNK.out.nc\", \"r\")'
        VAR1  ='VAR1 = addfile(\"' + job_path + '/ismip-hom-c/80km/data/ishom.c.80km.JFNK.out.nc\", \"r\")'
        png  = 'PNG = "' + ncl_path + '/ismipc80"'
        plot_ishomc = "ncl '" + VAR1 + "' '" + stock + "' '" + png + "' " + ishomc_plotfile

#TODO create an iteration plot and have that also in the html file 
        try:
                output = subprocess.call(plot_ishomc, shell=True)
        except:
                print "error creating ncl ismip hom c 80km plots"
                raise

# transferring ismipa pic to www file

        if (ncl_path + '/ismipc80.png'): 
                ishomcpic = "mv -f " + ncl_path + "/ismipc80.png" + " " + html_path + "/"
                try:
                        output = subprocess.call(ishomcpic, shell=True)
                except:
                        print "error moving ismip hom c 80km png file to www directory"
                        sys.exit(1)
                        raise

        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>ISMIP HOM C 80km </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark run</H4>\n')
        plot_file.write('<BR>\n')
        plot_file.write('<OBJECT data="ismipc80.png" type="image/png" width="1100" height="800" hspace=10 align=left alt="ISMIP HOM C Plots">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

        plot_file.write('</HTML>\n')
        plot_file.close()


