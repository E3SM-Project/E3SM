#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_utilities
import VV_outprocess
import VV_checks

# nomenclature for solver iteration values (1) d=dome, (2) d=diagnostic or e=evolving, (3) 30= size, 
# (4) 1,4,9,15 = processor count (5) b = bench if benchmark

def ddetails(solver_file,job_path,ncl_path,data_path,target_html):  # using data, fill the web page with info
        
        failedt_list = []
	solver_file.write('<HTML>\n')
	solver_file.write('<H3>Diagnostic Dome 30 Iteration Count Details:</H3>')
        solver_file.write('<H4>Eventually published in plot form</H4>')
	solver_file.write('<BR> \n')
# JFNK gnu 1 proc
#	print job_path + '/dome30/diagnostic/data/gnu.JFNK.1proc'

# Failure checking
        failedt1 = VV_checks.failcheck(job_path, '/dome30/diagnostic/data/gnu.JFNK.1proc')
        failedt_list.append(failedt1)

	solver_file.write('<H4>New Run: gnu.JFNK.1proc</H4>')
	procttl_dd301, nonlist_dd301, avg2_dd301, out_flag_dd301, ndd301_name, ldd301_name = VV_outprocess.jobprocess(job_path + '/dome30/diagnostic/data/gnu.JFNK.1proc', 'domed301')

	solver_file.write("Number of Processors = " + str(procttl_dd301[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
        VV_utilities.format(solver_file, nonlist_dd301)

# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR>\n')
	if out_flag_dd301 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_dd301)

# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run: gnu.JFNK.1proc</H4>')
	procttl_dd301b, nonlist_dd301b, avg2_dd301b, out_flag_dd301b, ndd301b_name, ldd301b_name = VV_outprocess.jobprocess(job_path + '/bench/dome30/diagnostic/data/gnu.JFNK.1proc', 'domed301b')

	solver_file.write("Number of Processors = " + str(procttl_dd301b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	VV_utilities.format(solver_file, nonlist_dd301b)
	solver_file.write('<BR>\n')
	if out_flag_dd301b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_dd301b)
	solver_file.write('<BR> \n')

# JFNK gnu 4 proc

# Failure checking
        failedt2 = VV_checks.failcheck(job_path, '/dome30/diagnostic/data/gnu.JFNK.4proc')
        failedt_list.append(failedt2)

        solver_file.write('<H4>New Run: gnu.JFNK.4proc</H4>')
	procttl_dd304, nonlist_dd304, avg2_dd304, out_flag_dd304, ndd304_name, ldd304_name = VV_outprocess.jobprocess(job_path + '/dome30/diagnostic/data/gnu.JFNK.4proc','domed304')

	solver_file.write("Number of Processors = " + str(procttl_dd304[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	VV_utilities.format(solver_file, nonlist_dd304)
	solver_file.write('<BR>\n')
	if out_flag_dd304 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_dd304)
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run: gnu.JFNK.4proc</H4>')
	procttl_dd304b, nonlist_dd304b, avg2_dd304b, out_flag_dd304b, ndd304b_name, ldd304b_name = VV_outprocess.jobprocess(job_path + '/bench/dome30/diagnostic/data/gnu.JFNK.4proc','domed304b')

	solver_file.write("Number of Processors = " + str(procttl_dd304b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	VV_utilities.format(solver_file, nonlist_dd304b)
	solver_file.write('<BR>\n')
	if out_flag_dd304b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_dd304b)
	solver_file.write('<BR> \n')

        if 1 in failedt_list:
            failedt = 1
        else:
            failedt = 0

        return failedt
def edetails(solver_file,job_path,ncl_path,data_path,target_html):  # using data, fill the web page with info

        failedt_list = []
	solver_file.write('<HTML>\n')
	solver_file.write('<H3>Evolving Dome 30 Iteration Count Details:</H3>')
        solver_file.write('<H4>Eventually published in plot form</H4>')
	solver_file.write('<BR> \n')
# JFNK gnu 9 proc
#	print job_path + '/dome30/evolving/data/gnu.JFNK.1proc'

# Failure checking
        failedt1 = VV_checks.failcheck(job_path, '/dome30/evolving/data/gnu.JFNK.9proc')
        failedt_list.append(failedt1)

        procttl_de309, nonlist_de309, avg2_de309, out_flag_de309, nde309_name, lde309_name = VV_outprocess.jobprocess(job_path + '/dome30/evolving/data/gnu.JFNK.9proc', 'domee309')
	procttl_de309b, nonlist_de309b, avg2_de309b, out_flag_de309b, nde309b_name, lde309b_name = VV_outprocess.jobprocess(job_path + '/bench/dome30/evolving/data/gnu.JFNK.9proc', 'domee309b')

# create iteration plots for proudction simulation

#	data_script=ncl_path + "/solver_dome30.ncl"

#	plot_dome30_data = "ncl 'nfile=\"" + data_path + "" + nde309_name + "\"'" + ' ' + \
#                     "'lfile=\"" + data_path + "" + lde309_name + "\"'" + ' ' + \
#                     "'PNG=\"" + ncl_path + "/dome30e_iter\"'" + ' ' + \
#                    data_script + ' ' + "1> /dev/null"

    
#	try:
#	       	output = subprocess.call(plot_dome30_data, shell=True)
#	except:
#        	print "error formatting iteration plot of dome30 evolving run"
#        	raise

#transferring iteration plot to www location

#	if (ncl_path + '/dome30e_iter.png'):
#        	iterpic = "mv -f " + ncl_path + "/dome30e_iter.png" + " " + target_html + "/"
#        	try:
#        	        output = subprocess.call(iterpic, shell=True)
#        	except:
#        	        print "error moving dome30 evolving iteration png file"
#        	        raise

# also display in list format
#        solver_file.write('<TABLE>\n')
#        solver_file.write('<TR>\n')
#        solver_file.write('<H4>Iteration Count for Evolving Dome30 Nonlinear and Linear Solver</H4>\n')
#        solver_file.write('<OBJECT data="dome30e_iter.png" type="image/png" width="1300" height="800" hspace=10 align=left alt="Solver Plots">\n')
#        solver_file.write('</OBJECT>\n')
#        solver_file.write('<TR>\n')
        solver_file.write('<BR>\n')
#        solver_file.write('</TABLE>\n')

	solver_file.write('<H4>New Run: gnu.JFNK.9proc</H4>')

	solver_file.write("Number of Processors = " + str(procttl_de309[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	VV_utilities.format(solver_file, nonlist_de309)
	solver_file.write('<BR>\n')
	if out_flag_de309 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_de309)
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run: gnu.JFNK.9proc</H4>')

	solver_file.write("Number of Processors = " + str(procttl_de309b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	VV_utilities.format(solver_file, nonlist_de309b)
	solver_file.write('<BR>\n')
	if out_flag_de309b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_de309b)
	solver_file.write('<BR> \n')

# JFNK gnu 15 proc
# Failure checking
        failedt2 = VV_checks.failcheck(job_path, '/dome30/evolving/data/gnu.JFNK.15proc')
        failedt_list.append(failedt2)

        solver_file.write('<H4>New Run: gnu.JFNK.15proc</H4>')
	procttl_de3015, nonlist_de3015, avg2_de3015, out_flag_de3015, nde3015_name, lde3015_name = VV_outprocess.jobprocess(job_path + '/dome30/evolving/data/gnu.JFNK.15proc', 'domee3015')

	solver_file.write("Number of Processors = " + str(procttl_de3015[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	VV_utilities.format(solver_file, nonlist_de3015)
	solver_file.write('<BR>\n')
	if out_flag_de3015 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_de3015)
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run: gnu.JFNK.15proc</H4>')
	procttl_de3015b, nonlist_de3015b, avg2_de3015b, out_flag_de3015b, nde3015b_name, lde3015b_name = VV_outprocess.jobprocess(job_path + '/bench/dome30/evolving/data/gnu.JFNK.15proc', 'domee3015b')

	solver_file.write("Number of Processors = " + str(procttl_de3015b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	VV_utilities.format(solver_file, nonlist_de3015b)
	solver_file.write('<BR>\n')
	if out_flag_de3015b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	VV_utilities.format(solver_file, avg2_de3015b)
	solver_file.write('<BR> <BR>\n')

	solver_file.write('</HTML>\n')
	solver_file.close()

        if 1 in failedt_list:
            failedt = 1
        else:
            failedt = 0

        return failedt

def dplot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info
        
        noplot = 0

        tmpath = job_path + '/dome30/diagnostic/data/dome.1.nc'
        tmpath2 = job_path + '/dome30/diagnostic/data/dome.4.nc'
        if VV_utilities.emptycheck(tmpath) == 0 or VV_utilities.emptycheck(tmpath2) == 0:
                noplot = 1

	plot_file.write('<HTML>\n')
	plot_file.write('<H3>Diagnostic Dome 30 Plot Details:</H3>')
	
        dome30_plotfile=''+ ncl_path + '/dome30d.ncl'
	stock='STOCK = addfile(\"'+ job_path + '/bench/dome30/diagnostic/data/dome.out.nc\", \"r\")'
	VAR1  ='VAR1 = addfile(\"' + job_path + '/dome30/diagnostic/data/dome.1.nc\", \"r\")'
	VAR4  ='VAR4 = addfile(\"' + job_path + '/dome30/diagnostic/data/dome.4.nc\", \"r\")'
	png  = 'PNG = "' + ncl_path + '/dome30d"'
        plot_dome30d = "ncl '" + VAR1 + "' '" + VAR4 + \
                           "' '" + stock + "' '" + png + "' " + dome30_plotfile 
#        print plot_dome30

        try:
                output = subprocess.call(plot_dome30d, shell=True)
        except:
                print "error creating ncl diagnostic dome30 plots"
                raise

# transferring dome30 pic to www file

        if (ncl_path + '/dome30d.png'):
                dome30dpic = "mv -f " + ncl_path + "/dome30d.png" + " " + html_path + "/"
        	try:
                        output = subprocess.call(dome30dpic, shell=True)
        	except:
                	print "error moving diagnostic dome30 png file to www directory"
                        sys.exit(1)
                	raise

        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>Diagnostic Dome 30 </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark for 1 and 4 processors, Velocity Norm </H4>\n')
        plot_file.write('<OBJECT data="dome30d.png" type="image/png" width="1100" height="800" hspace=5 align=left alt="Dome 30 Plots">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

	plot_file.write('</HTML>\n')
	plot_file.close()

def eplot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

	plot_file.write('<HTML>\n')
	plot_file.write('<H3>Evolving Dome 30 Plot Details:</H3>')
	dome30_plotfile=''+ ncl_path + '/dome30e.ncl'
	stock='STOCK = addfile(\"'+ job_path + '/bench/dome30/evolving/data/dome.out.nc\", \"r\")'
	VAR9  ='VAR9 = addfile(\"' + job_path + '/dome30/evolving/data/dome.9.nc\", \"r\")'
	VAR15  ='VAR15 = addfile(\"' + job_path + '/dome30/evolving/data/dome.15.nc\", \"r\")'
	png  = 'PNG = "' + ncl_path + '/dome30e"'
        plot_dome30e = "ncl '" + VAR9 + "' '" + VAR15 + \
                           "' '" + stock + "' '" + png + "' " + dome30_plotfile 
#        print plot_dome30

        try:
                output = subprocess.call(plot_dome30e, shell=True)
        except:
                print "error creating ncl evolving dome30 velocity and thickness plot"
                raise

# transferring dome30 pic to www file

        if (ncl_path + '/dome30e.png'):
        	dome30epic = "mv -f " + ncl_path + "/dome30e.png" + " " + html_path + "/"
        	try:
                	output = subprocess.call(dome30epic, shell=True)
        	except:
                	print "error moving evolving dome30 vel/thk png file to www directory"
                        sys.exit(1)
                	raise

        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>Evolving Dome 30 </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark for 9 and 15 processors, Vel Norm and Thickness </H4>\n')
        plot_file.write('<OBJECT data="dome30e.png" type="image/png" width="1100" height="800" hspace=5 align=left alt="Evolving Dome 30 Velocity and Thickness Plots">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

	dome30_plotfile=''+ ncl_path + '/dome30t.ncl'
	png  = 'PNG = "' + ncl_path + '/dome30t"'
        plot_dome30 = "ncl '" + VAR9 + "' '" + VAR15 + \
                           "' '" + stock + "' '" + png + "' " + dome30_plotfile 
        try:
                output = subprocess.call(plot_dome30, shell=True)
        except:
                print "error creating ncl evolving dome30 temperature plots"
                raise

# transferring dome30 pic to www file

        if (ncl_path + '/dome30t.png'):
        	dome30pic = "mv -f " + ncl_path + "/dome30t.png" + " " + html_path + "/"
        	try:
                	output = subprocess.call(dome30pic, shell=True)
        	except:
                	print "error moving evolving dome30 temp png file to www directory"
                        sys.exit(1)
                	raise

        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark for 9 and 15 processors, Temperature </H4>\n')
        plot_file.write('<OBJECT data="dome30t.png" type="image/png" width="1100" height="800" hspace=5 align=left alt="Evolving Dome 30 Temperature Plots">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

	plot_file.write('</HTML>\n')
	plot_file.close()

