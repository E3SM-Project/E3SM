#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_outprocess

# nomenclature for solver iteration values (1) d=dome, (2) d=diagnostic or e=evolving, (3) 30= size, 
# (4) 1,4,9,15 = processor count (5) b = bench if benchmark

def ddetails(solver_file,job_path,ncl_path,data_path,target_html):  # using data, fill the web page with info

	solver_file.write('<HTML>\n')
	solver_file.write('<H3>Diagnostic Dome 30 Iteration Count Details:</H3>')

	solver_file.write('<BR> \n')
# JFNK gnu 1 proc
#	print job_path + '/dome30/diagnostic/data/gnu.JFNK.1proc'

	solver_file.write('<H4>New Run:</H4>')
	procttl_dd301, nonlist_dd301,avg2_dd301,out_flag_dd301,ndd301_name,ldd301_name = VV_outprocess.jobprocess(job_path + '/dome30/diagnostic/data/gnu.JFNK.1proc', 'domed301')

	solver_file.write("Number of Processors = " + str(procttl_dd301[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_dd301:
		solver_file.write(str(item) + ", ")
# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR>\n')
	if out_flag_dd301 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_dd301:
		solver_file.write(str(item) + ", ")
# print this to an array that ncl can use for plotting
#                print item
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_dd301b, nonlist_dd301b,avg2_dd301b,out_flag_dd301b,ndd301b_name,ldd301b_name = VV_outprocess.jobprocess(job_path + '/dome30/diagnostic/bench/gnu.JFNK.1proc', 'domed301b')

	solver_file.write("Number of Processors = " + str(procttl_dd301b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_dd301b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_dd301b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_dd301b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

# JFNK gnu 4 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_dd304, nonlist_dd304,avg2_dd304,out_flag_dd304,ndd304_name,ldd304_name = VV_outprocess.jobprocess(job_path + '/dome30/diagnostic/data/gnu.JFNK.4proc','domed304')

	solver_file.write("Number of Processors = " + str(procttl_dd304[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_dd304:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_dd304 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_dd304:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_dd304b, nonlist_dd304b,avg2_dd304b,out_flag_dd304b,ndd304b_name,ldd304b_name = VV_outprocess.jobprocess(job_path + '/dome30/diagnostic/bench/gnu.JFNK.4proc','domed304b')

	solver_file.write("Number of Processors = " + str(procttl_dd304b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_dd304b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_dd304b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) THAT FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_dd304b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

def edetails(solver_file,job_path,ncl_path,data_path,target_html):  # using data, fill the web page with info

	solver_file.write('<HTML>\n')
	solver_file.write('<H3>Evolving Dome 30 Iteration Count Details:</H3>')

	solver_file.write('<BR> \n')
# JFNK gnu 9 proc
#	print job_path + '/dome30/evolving/data/gnu.JFNK.1proc'
	procttl_de309, nonlist_de309,avg2_de309,out_flag_de309,nde309_name,lde309_name = VV_outprocess.jobprocess(job_path + '/dome30/evolving/data/gnu.JFNK.9proc', 'domee309')
	procttl_de309b, nonlist_de309b,avg2_de309b,out_flag_de309b,nde309b_name,lde309b_name = VV_outprocess.jobprocess(job_path + '/dome30/evolving/bench/gnu.JFNK.9proc', 'domee309b')

# create iteration plots for proudction simulation

	data_script=ncl_path + "/solver_dome30.ncl"

	plot_dome30_data = "ncl 'nfile=\"" + data_path + "" + nde309_name + "\"'" + ' ' + \
                     "'lfile=\"" + data_path + "" + lde309_name + "\"'" + ' ' + \
                     "'PNG=\"" + ncl_path + "/dome30e_iter\"'" + ' ' + \
                    data_script + ' ' + "1> /dev/null"

	try:
	       	output = subprocess.call(plot_dome30_data, shell=True)
	except:
        	print "error formatting iteration plot of dome30 evolving run"
        	raise

#transferring iteration plot to www location

	if (ncl_path + '/dome30e_iter.png'):
        	iterpic = "mv -f " + ncl_path + "/dome30e_iter.png" + " " + target_html + "/"
        	try:
        	        output = subprocess.call(iterpic, shell=True)
        	except:
        	        print "error moving dome30 evolving iteration png file"
        	        raise

        solver_file.write('<TABLE>\n')
        solver_file.write('<TR>\n')
        solver_file.write('<H4>Iteration Count for Evolving Dome30 Nonlinear and Linear Solver</H4>\n')
        solver_file.write('<OBJECT data="dome30e_iter.png" type="image/png" width="1300" height="800" hspace=10 align=left alt="Solver Plots">\n')
        solver_file.write('</OBJECT>\n')
        solver_file.write('<TR>\n')
        solver_file.write('<BR>\n')
        solver_file.write('</TABLE>\n')

# JFNK gnu 9 proc
	solver_file.write('<H4>New Run:</H4>')

	solver_file.write("Number of Processors = " + str(procttl_de309[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_de309:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_de309 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_de309:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')

	solver_file.write("Number of Processors = " + str(procttl_de309b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_de309b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_de309b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_de309b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

# JFNK gnu 15 proc
	solver_file.write('<H4>New Run:</H4>')
	procttl_de3015, nonlist_de3015,avg2_de3015,out_flag_de3015,nde3015_name,lde3015_name = VV_outprocess.jobprocess(job_path + '/dome30/evolving/data/gnu.JFNK.15proc', 'domee3015')

	solver_file.write("Number of Processors = " + str(procttl_de3015[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_de3015:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_de3015 == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_de3015:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> \n')

	solver_file.write('<H4>Benchmark Run:</H4>')
	procttl_de3015b, nonlist_de3015b,avg2_de3015b,out_flag_de3015b,nde3015b_name,lde3015b_name = VV_outprocess.jobprocess(job_path + '/dome30/evolving/bench/gnu.JFNK.15proc', 'domee3015b')

	solver_file.write("Number of Processors = " + str(procttl_de3015b[-1]) + "<BR>\n")
	solver_file.write("Number of Nonlinear Iterations = ")
	for item in nonlist_de3015b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR>\n')
	if out_flag_de3015b == 1:
		solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
	solver_file.write("Average Number of Linear Iterations per Time-Step = ")
	for item in avg2_de3015b:
		solver_file.write(str(item) + ", ")
	solver_file.write('<BR> <BR>\n')

	solver_file.write('</HTML>\n')
	solver_file.close()

def dplot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

	plot_file.write('<HTML>\n')
	plot_file.write('<H3>Diagnostic Dome 30 Plot Details:</H3>')
	dome30_plotfile=''+ ncl_path + '/dome30d.ncl'
	stock='STOCK = addfile(\"'+ job_path + '/dome30/diagnostic/bench/dome.out.nc\", \"r\")'
	VAR1  ='VAR1 = addfile(\"' + job_path + '/dome30/diagnostic/data/dome.1.nc\", \"r\")'
	VAR4  ='VAR4 = addfile(\"' + job_path + '/dome30/diagnostic/data/dome.4.nc\", \"r\")'
	png  = 'PNG = "' + ncl_path + '/dome30d"'
        plot_dome30 = "ncl '" + VAR1 + "' '" + VAR4 + \
                           "' '" + stock + "' '" + png + "' " + dome30_plotfile 
#        print plot_dome30

        try:
                output = subprocess.call(plot_dome30, shell=True)
        except:
                print "error creating ncl diagnostic dome30 plots"
                raise

# transferring dome30 pic to www file

        if (ncl_path + '/dome30d.png'):
        	dome30pic = "mv -f " + ncl_path + "/dome30d.png" + " " + html_path + "/"
        	try:
                	output = subprocess.call(dome30pic, shell=True)
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
	stock='STOCK = addfile(\"'+ job_path + '/dome30/evolving/bench/dome.out.nc\", \"r\")'
	VAR9  ='VAR9 = addfile(\"' + job_path + '/dome30/evolving/data/dome.9.nc\", \"r\")'
	VAR15  ='VAR15 = addfile(\"' + job_path + '/dome30/evolving/data/dome.15.nc\", \"r\")'
	png  = 'PNG = "' + ncl_path + '/dome30e"'
        plot_dome30 = "ncl '" + VAR9 + "' '" + VAR15 + \
                           "' '" + stock + "' '" + png + "' " + dome30_plotfile 
#        print plot_dome30

        try:
                output = subprocess.call(plot_dome30, shell=True)
        except:
                print "error creating ncl evolving dome30 velocity and thickness plot"
                raise

# transferring dome30 pic to www file

        if (ncl_path + '/dome30e.png'):
        	dome30pic = "mv -f " + ncl_path + "/dome30e.png" + " " + html_path + "/"
        	try:
                	output = subprocess.call(dome30pic, shell=True)
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

def case(plot_file,confpath,benchpath,ncl_path,html_path):  # using data, fill the web page with info

# grab and print config file information
# this function allows creation of nested dictionaries on the fly (like PERL autovivification)
	def makehashdome():
		return collections.defaultdict(makehashdome)
	datadome = makehashdome()

	def makehashbench():
		return collections.defaultdict(makehashbench)
	databench = makehashbench()

	con_flag = False
	keywords = ('parameters', 'CF output', 'grid', 'time', 'options', 'ho_options')
	variable, value = '', ''
	for kw in keywords:
	       	try:
       		 	configlog = open(confpath, 'r')
		
		except:
       	         	print "error reading" + configlog
                	sys.exit(1)
                	raise
	       	try:
       		 	benchlog = open(benchpath, 'r')
		
		except:
       	         	print "error reading" + benchlog 
                	sys.exit(1)
                	raise

                for cfline in configlog:
                        if cfline.startswith('[' + kw + ']'):
                               con_flag = True
                               continue
                        if cfline.startswith('['):
                                con_flag = False
                                continue
                        if con_flag == True:
                                line = cfline.strip('\r\n')
                                if '#' in line:
                                        tmp = line.split('#')
                                        line = tmp[0]
                                if line == '':
                                        continue
                                if line.endswith('='):
                                        variable, junk = line.split()
                                        value = ''
                                else:
                                        variable, value = line.split('=')
                                variable = variable.strip()
                                value = value.strip()
                                datadome[kw][variable] = value

                configlog.close()

                for cfline in benchlog:
                        if cfline.startswith('[' + kw + ']'):
                               con_flag = True
                               continue
                        if cfline.startswith('['):
                                con_flag = False
                                continue
                        if con_flag == True:
                                line = cfline.strip('\r\n')
                                if '#' in line:
                                        tmp = line.split('#')
                                        line = tmp[0]
                                if line == '':
                                        continue
                                if line.endswith('='):
                                        variable, junk = line.split()
                                        value = ''
                                else:
                                        variable, value = line.split('=')
                                variable = variable.strip()
                                value = value.strip()
                                databench[kw][variable] = value

                configlog.close()

#Calculate number of time steps: dome
	if datadome['time']['tend'] and datadome['time']['tstart']:
        	diff = float(datadome['time']['tend']) - float(datadome['time']['tstart'])
        	timestp = diff / float(datadome['time']['dt'])
#Calculate number of time steps: bench
	if databench['time']['tend'] and databench['time']['tstart']:
        	diff = float(databench['time']['tend']) - float(databench['time']['tstart'])
        	timestp_bench = diff / float(databench['time']['dt'])

#Put settings on website, check if they match the benchmark settings
	plot_file.write('<HTML>\n')
	plot_file.write('<H3>Case Details:</H3>')
        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>Configure Settings </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')

        plot_file.write('<H4>Configure file Settings </H4>\n')
        if datadome['CF output']['variables']:
                plot_file.write('Output available from test run: ' + datadome['CF output']['variables'] + "<BR>\n")

#Grid Size/spacing
        plot_file.write("Grid Size (vert by ew by ns): " + datadome['grid']['upn'] + "x" + datadome['grid']['ewn'] + "x" + datadome['grid']['nsn'] + "<BR>\n")
        if datadome['grid']['upn'] != databench['grid']['upn']:
        	plot_file.write('<FONT COLOR="red">vertical dimension (upn) is different than benchmark run: ' + datadome['grid']['upn'] + ' versus ' + databench['grid']['upn'] + '</FONT><BR>\n')
        if datadome['grid']['ewn'] != databench['grid']['ewn']:
        	plot_file.write('<FONT COLOR="red">east/west dimension (ewn) is different than benchmark run: ' + datadome['grid']['ewn'] + ' versus ' + databench['grid']['ewn'] + '</FONT><BR>\n')
        if datadome['grid']['nsn'] != databench['grid']['nsn']:
        	plot_file.write('<FONT COLOR="red">north/south dimension (nsn) is different than benchmark run: ' + datadome['grid']['nsn'] + ' versus ' + databench['grid']['nsn'] + '</FONT><BR>\n')

        plot_file.write("Grid Spacing (ew by ns): " + datadome['grid']['dew'] + "x" + datadome['grid']['dns'] + "<BR>\n")
        if datadome['grid']['dew'] != databench['grid']['dew']:
        	plot_file.write('<FONT COLOR="red">east/west spacing (dew) is different than benchmark run: ' + datadome['grid']['dew'] + ' versus ' + databench['grid']['dew'] + '</FONT><BR>\n')
        if datadome['grid']['dns'] != databench['grid']['dns']:
        	plot_file.write('<FONT COLOR="red">north/south spacing (dns) is different than benchmark run: ' + datadome['grid']['dns'] + ' versus ' + databench['grid']['dns'] + '</FONT><BR>\n')

#Simulation length and time steps
        plot_file.write("Start/End Time: " + datadome['time']['tstart'] + "," + datadome['time']['tend'] + ", Number of time steps = " + str(timestp) + "<BR>\n")
        if datadome['time']['tstart'] != databench['time']['tstart']:
        	plot_file.write('<FONT COLOR="red">Start time is different than benchmark run: ' + datadome['time']['tstart'] + ' versus ' + databench['time']['tstart'] + '</FONT><BR>\n')
        if datadome['time']['tend'] != databench['time']['tend']:
        	plot_file.write('<FONT COLOR="red">Start time is different than benchmark run: ' + datadome['time']['tend'] + ' versus ' + databench['time']['tend'] + '</FONT><BR>\n')
        if str(timestp) != str(timestp_bench):
        	plot_file.write('<FONT COLOR="red">Number of timesteps is different than benchmark run: ' + str(timestp) + ' versus ' + str(timestp_bench) + '</FONT><BR>\n')

#Parameter settings
        plot_file.write('<BR\n>')
        plot_file.write("Parameters<BR>\n")
        plot_file.write("------------------------<BR>\n")
        if datadome['parameters']['flow_factor']:
        	if datadome['parameters']['flow_factor'] == databench['parameters']['flow_factor']:
                	plot_file.write('flow_factor = ' + datadome['parameters']['flow_factor'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">flow_factor = ' + datadome['parameters']['flow_factor'] + ' different than benchmark value:  ' + databench['parameters']['flow_factor'] + '</FONT><BR>\n')
        if datadome['parameters']['ice_limit']:
        	if datadome['parameters']['ice_limit'] == databench['parameters']['ice_limit']:
                	plot_file.write('ice_limit = ' + datadome['parameters']['ice_limit'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">Ice Limit = ' + datadome['parameters']['ice_limit'] + ' different than benchmark value:  ' + databench['parameters']['ice_limit'] + '</FONT><BR>\n')
        if datadome['parameters']['default_flwa']:
        	if datadome['parameters']['default_flwa'] == databench['parameters']['default_flwa']:
                	plot_file.write('default_flwa = ' + datadome['parameters']['default_flwa'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">default_flwa = ' + datadome['parameters']['default_flwa'] + ' different than benchmark value:  ' + databench['parameters']['default_flwa'] + '</FONT><BR>\n')

#Options settings
        plot_file.write('<BR\n>')
        plot_file.write("Options<BR>\n")
        plot_file.write("------------------------<BR>\n")
        if datadome['options']['dycore']:
        	if datadome['options']['dycore'] == databench['options']['dycore']:
                	plot_file.write('dycore = ' + datadome['options']['dycore'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">dycore = ' + datadome['options']['dycore'] + ' different than benchmark value:  ' + databench['options']['dycore'] + '</FONT><BR>\n')
        if datadome['options']['flow_law']:
        	if datadome['options']['flow_law'] == databench['options']['flow_law']:
                	plot_file.write('flow_law = ' + datadome['options']['flow_law'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">flow_law = ' + datadome['options']['flow_law'] + ' different than benchmark value:  ' + databench['options']['flow_law'] + '</FONT><BR>\n')
        if datadome['options']['evolution']:
        	if datadome['options']['evolution'] == databench['options']['evolution']:
                	plot_file.write('evolution = ' + datadome['options']['evolution'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">evolution = ' + datadome['options']['evolution'] + ' different than benchmark value:  ' + databench['options']['evolution'] + '</FONT><BR>\n')
        if datadome['options']['temperature']:
        	if datadome['options']['temperature'] == databench['options']['temperature']:
                	plot_file.write('temperature = ' + datadome['options']['temperature'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">temperature = ' + datadome['options']['temperature'] + ' different than benchmark value:  ' + databench['options']['temperature'] + '</FONT><BR>\n')

#HO Options settings
        plot_file.write('<BR\n>')
        plot_file.write("HO Options<BR>\n")
        plot_file.write("------------------------<BR>\n")
        if datadome['ho_options']['diagnostic_scheme']:
        	if datadome['ho_options']['diagnostic_scheme'] == databench['ho_options']['diagnostic_scheme']:
                	plot_file.write('diagnostic_scheme = ' + datadome['ho_options']['diagnostic_scheme'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">diagnostic_scheme = ' + datadome['ho_options']['diagnostic_scheme'] + ' different than benchmark value:  ' + databench['ho_options']['diagnostic_scheme'] + '</FONT><BR>\n')
        if datadome['ho_options']['which_ho_babc']:
        	if datadome['ho_options']['which_ho_babc'] == databench['ho_options']['which_ho_babc']:
                	plot_file.write('which_ho_babc = ' + datadome['ho_options']['which_ho_babc'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">which_ho_babc = ' + datadome['ho_options']['which_ho_babc'] + ' different than benchmark value:  ' + databench['ho_options']['which_ho_babc'] + '</FONT><BR>\n')
        if datadome['ho_options']['which_ho_efvs']:
        	if datadome['ho_options']['which_ho_efvs'] == databench['ho_options']['which_ho_efvs']:
                	plot_file.write('which_ho_efvs = ' + datadome['ho_options']['which_ho_efvs'] + "<BR>\n")
		else:
        		plot_file.write('<FONT COLOR="red">which_ho_efvs = ' + datadome['ho_options']['which_ho_efvs'] + ' different than benchmark value:  ' + databench['ho_options']['which_ho_efvs'] + '</FONT><BR>\n')
#        if datadome['ho_options']['which_ho_nonlinear']:
#        	if datadome['ho_options']['which_ho_nonlinear'] == databench['ho_options']['which_ho_nonlinear']:
#                	plot_file.write('which_ho_nonlinear = ' + datadome['ho_options']['which_ho_nonlinear'] + "<BR>\n")
#		else:
#        		plot_file.write('<FONT COLOR="red">which_ho_nonlinear = ' + datadome['ho_options']['which_ho_nonlinear'] + ' different than benchmark value:  ' + databench['ho_options']['which_ho_nonlinear'] + '</FONT><BR>\n')
        plot_file.write('<BR\n>')

	plot_file.write('</HTML>\n')
	plot_file.close()
