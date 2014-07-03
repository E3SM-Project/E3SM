
import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_outprocess
import VV_testsuite

user = os.environ['USER']

#enable use input for filepaths from the command line
usage_string = "%prog [options]"
parser = OptionParser(usage=usage_string)
parser.add_option('-p', '--directory', action='store', type='string', dest='directory_path', \
                  metavar='PATH', help='path where this directory is located')
parser.add_option('-j', '--html', action='store', type='string', dest='html_path', \
                  metavar='PATH', help='path where the html directory is located')
parser.add_option('-l', '--link', action='store', type='string', dest='html_link', \
                  metavar='PATH', help='location of website for viewing set by user')
parser.add_option('-k', '--ncl', action='store', type='string', dest='ncl_path', \
                  metavar='PATH', help='path where the ncl directory is located')
parser.add_option('-d', '--data', action='store', type='string', dest='data_path', \
                  metavar='PATH', help='path where the solver data directory is located')
parser.add_option('-c', '--config', action='store', type='string', dest='config_file', \
                  metavar='FILE', help='the config file python will parse through')
parser.add_option('-t', '--test', action='store', type='string', dest='test_suite', \
                  metavar='TEST', help='path to location of test suite')
parser.add_option('-o', '--output', action='store', type='string', dest='output_file', \
                  metavar='FILE', help='the job output file python will parse through')
parser.add_option('-s', '--stocknc', action='store', type='string', dest='stock_netcdf_file', \
                  metavar='FILE', help='the stock NETCDF file that the ncl script will read')
parser.add_option('-n', '--variablenc', action='store', type='string', dest='variable_netcdf_file', \
                  metavar='FILE', help='the variable NETCDF file that the ncl script will read')
parser.add_option('-i', '--timestamp', action='store', type='string', dest='time_stamp', \
                  metavar='FILE', help='the current time to record in the web output')
parser.add_option('-m', '--comment', action='store', type='string', dest='comment', \
                  metavar='FILE', help='information about the test case for user reference')
parser.add_option('-u', '--username', action='store', type='string', dest='username', \
                  metavar='FILE', help='username used to create subdirectory of web pages of output')
parser.add_option('-g', '--gis_prod', action='store_true', dest='gis_prod', \
                  help='include flag to run the GIS production analysis')
#parser.add_option('-a', '--ant_prod', action='store_true', dest='ant_prod', \
#                  help='include flag to run the ANT production analysis')

#parse the command line options and arguments and store in lists
(options, args) = parser.parse_args()

if (options.gis_prod):
	if (options.stock_netcdf_file):
		stock_nc = options.stock_netcdf_file 
	else:
		print "need a benchmark GIS file for production analysis"
	        exit()

	if (options.variable_netcdf_file):
        	variable_nc = options.variable_netcdf_file
	else:
		print "need a production GIS file for analysis"
       		exit()

	if (options.output_file):
		gis_output = options.output_file 
	else:
		print "no GIS output file provided, so no solver statistics will be provided"

else:
        print "not performing GIS production analysis"

if (options.time_stamp):
	time_stamp = options.time_stamp 
else:
	print "no time given for website"
	time_stamp = " "

if (options.comment):
	comment = options.comment
else:
	print "no comments about test case given"
	comment = " "

if (options.username):

	print 'placing HTML files in the ' + options.username + ' subdirectory (check permissions)'
 	target_html = options.html_path + '/' + options.username

else:

	print 'no username specified, placing HTML files in the main html directory'
	target_html = options.html_path 

#remove html files previously used in specified subdirectory

try:
	os.remove(target_html + '/GIS-main-diag.html')
except OSError as o:
      	if o.errno == 2:
		print "recreating GIS-main-diag.html in " + options.username + " subdirectory"
	else:
		raise
try:
       	os.remove(target_html + '/test_suite.html')
except OSError as o:
	if o.errno == 2:
		print "recreating test suite in " + options.username + " subdirectory"
	else:
		raise
try:
       	os.remove(target_html + '/GIS-con-diag.html')
except OSError as o:
	if o.errno == 2:
		print "recreating GIS-con-diag.html in " + options.username + " subdirectory"
	else:
		raise
try:
       	os.remove(target_html + '/GIS-out-diag.html')
except OSError as o:
	if o.errno == 2:
		print "recreating GIS-out-diag.html in " + options.username + " subdirectory"
	else:
		raise
try:
       	os.remove(target_html + '/GIS-plot-diag.html')
except OSError as o:
	if o.errno == 2:
		print "recreating GIS-plot-diag.html in " + options.username + " subdirectory"
	else:
		raise

#read the configure file of the production run to provide simulation details

# this function allows creation of nested dictionaries on the fly (like PERL autovivification)
def makehash():
        return collections.defaultdict(makehash)
data = makehash()

con_flag = False
keywords = ('parameters', 'CF output', 'grid', 'time', 'options')
variable, value = '', ''

#open and grab info from configure file 
#TODO can grab input and output file information from this file, rather than having user specify in bash script
if options.gis_prod and options.config_file:
	for kw in keywords:
		try:
			configlog = open(options.config_file, "r")
		except:
			print "error reading configure file, or no file specified"
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
				data[kw][variable] = value

		configlog.close()
	
#Calculate number of time steps
	if data['time']['tend'] and data['time']['tstart']:
		diff = float(data['time']['tend']) - float(data['time']['tstart'])
		timestp = diff / float(data['time']['dt'])

# create production plots for analysis
if options.gis_prod:

	ncl_script=options.ncl_path + "/tri_ncl_script.ncl"

	plot_gis_thk = "ncl 'VAR=addfile(\"" + variable_nc + "\", \"r\")'" + ' ' + \
		    "'STOCK=addfile(\"" + stock_nc + "\",\"r\")'" + ' ' + \
		    "'PNG=\"" + options.ncl_path + "/gis5km_thk\"'" + ' ' + \
                    ncl_script + ' ' + "1> /dev/null"

#print plot_gis_thk

	try:
		output = subprocess.call(plot_gis_thk, shell=True)
	except:
		print "error formatting thickness plot of production run"
		raise

#transferring thickness plot to www location

	if (options.ncl_path + '/gis5km_thk.png'):
        	thkpic = "mv -f " + options.ncl_path + "/gis5km_thk.png" + " " + target_html + "/"
        	try:
                	output = subprocess.call(thkpic, shell=True)
        	except:
                	print "error moving thk png file"
                	raise

#transferring cover picture to www file

if (options.ncl_path + '/alaska_pic.png'):
        alaskapic = "cp -f " + options.ncl_path + "/alaska_pic.png" + " " + target_html + "/"
        try:
                output = subprocess.call(alaskapic, shell=True)
        except:
                print "error moving cover picture"
                raise

#writing the main HTML page

file = open(target_html + '/GIS-main-diag.html', 'w')

file.write('<HTML><HEAD>\n')
file.write('<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">')
file.write('<title>Land Ice Verification and Validation toolkit</title>')
file.write('<link href="style.css" rel="stylesheet" type="text/css">')
file.write('<BODY>\n')
file.write(' <div id="container"> <div id="content">')
#file.write('<BODY bgcolor="white">\n')
file.write('<P>\n')
file.write('<OBJECT data="alaska_pic.png" type="image/png" width="400" height="300" hspace=10 align=left alt="ice sheet pic">\n')
file.write('</OBJECT>\n')
file.write('<FONT color=blue><B>\n')
file.write('Land Ice Validation package  </B></FONT> <BR>\n')
file.write('Performed on ' + options.time_stamp + '<BR>\n')
file.write('Test case run by: ' + options.username + '<BR>\n')
file.write('Details:' + options.comment + '<BR>\n')
file.write('</P>\n')
file.write('<BR clear=left>\n')
file.write('<BR>\n')
file.write('<HR noshade size=2 size="100%">\n')
file.write('<TH ALIGN=LEFT><A HREF="test_suite.html">Basic Test Suite Diagnostics</A>\n')
file.write('<BR>\n')
file.write('<BR>\n')

if options.gis_prod:
	file.write('<TH ALIGN=LEFT><A HREF="GIS-con-diag.html">Production Configure Diagnostics</A>\n')
	file.write('<BR>\n')
	file.write('<BR>\n')
	file.write('<TH ALIGN=LEFT><A HREF="GIS-out-diag.html">Production Output Diagnostics</A>\n')
	file.write('<BR>\n')
	file.write('<BR>\n')
	file.write('<TH ALIGN=LEFT><A HREF="GIS-plot-diag.html">Ice Thickness</A>\n')
	file.write('<BR>\n <BR>\n <BR> \n')

file.write('<h4> For Additional Information: </h4> <p>')
file.write(' Kate Evans <br>')
file.write('Oak Ridge National Laboratory<br>')
file.write('1 Bethel Valley Road <br>')
file.write('Oak Ridge, Tennessee 37831-6015 <br>')
file.write('Email: evanskj at ornl dot gov <br> </p>')

#create all the test suite diagnostics pages 

if options.test_suite:

	test_file = open(target_html + '/test_suite.html', 'w')
	descript_file = open(target_html + '/test_descript.html', 'w')
# diagnostic dome case
	dome30d_file = open(target_html + '/dome30d_details.html', 'w')
	dome30d_case = open(target_html + '/dome30d_case.html', 'w')
	dome30d_plot = open(target_html + '/dome30d_plot.html', 'w')
# evolving dome case
	dome30e_file = open(target_html + '/dome30e_details.html', 'w')
	dome30e_case = open(target_html + '/dome30e_case.html', 'w')
	dome30e_plot = open(target_html + '/dome30e_plot.html', 'w')
# circular shelf case
	circ_file = open(target_html + '/circ_details.html', 'w')
	circ_case = open(target_html + '/circ_case.html', 'w')
	circ_plot = open(target_html + '/circ_plot.html', 'w')
# confined shelf case
	conf_file = open(target_html + '/conf_details.html', 'w')
	conf_case = open(target_html + '/conf_case.html', 'w')
	conf_plot = open(target_html + '/conf_plot.html', 'w')
# ismip hom a 80km case
	ishoma80_file = open(target_html + '/ishoma80_details.html', 'w')
	ishoma80_case = open(target_html + '/ishoma80_case.html', 'w')
	ishoma80_plot = open(target_html + '/ishoma80_plot.html', 'w')
# ismip hom c 80km case
	ishomc80_file = open(target_html + '/ishomc80_details.html', 'w')
	ishomc80_case = open(target_html + '/ishomc80_case.html', 'w')
	ishomc80_plot = open(target_html + '/ishomc80_plot.html', 'w')
# 10km GIS case
	gis10_file = open(target_html + '/gis10_details.html', 'w')
	gis10_case = open(target_html + '/gis10_case.html', 'w')
	gis10_plot = open(target_html + '/gis10_plot.html', 'w')

# TODO create a list of the html files of the included cases, then pass through the testsuite.web call

#path to python code to create all the test suite pages and data
	reg_test = options.test_suite + "/reg_test"

	VV_testsuite.web(descript_file,test_file,dome30d_file,dome30d_case,dome30d_plot, \
		dome30e_file,dome30e_case,dome30e_plot, \
		circ_file,circ_case,circ_plot,conf_file,conf_case,conf_plot, \
		ishoma80_file,ishoma80_case,ishoma80_plot,ishomc80_file,ishomc80_case,ishomc80_plot,
		gis10_file,gis10_case,gis10_plot, \
		reg_test,options.ncl_path,options.data_path,target_html)

#create www page with config information

if options.gis_prod:
	con_file = open(target_html + '/GIS-con-diag.html', 'w')

	con_file.write('<HTML>\n')
	con_file.write('<TITLE>GIS Configure Diagnostics</TITLE>\n')
	con_file.write('<H2>Configure File Diagnostics</H2>')
	con_file.write("CF Output<BR>\n")
	con_file.write("-----------------------<BR>\n")
	if data['CF output']['variables']:
		con_file.write('variables = ' + data['CF output']['variables'] + "<BR>\n")
	con_file.write('<BR\n>')
	con_file.write("Grid<BR>\n")
	con_file.write("-----------------------<BR>\n")
	if data['grid']['upn']:	
		con_file.write("upn = " + data['grid']['upn'] + "<BR>\n")
	if data['grid']['ewn']:	
		con_file.write("ewn = " +  data['grid']['ewn'] + "<BR>\n")
	if data['grid']['nsn']:	
		con_file.write("nsn = " + data['grid']['nsn'] + "<BR>\n")
	if data['grid']['dew']:	
		con_file.write("dew = " + data['grid']['dew'] + "<BR>\n")
	if data['grid']['dns']:	
		con_file.write("dns = " + data['grid']['dns'] + "<BR>\n")
	con_file.write('<BR\n>')
	con_file.write("Time<BR>\n")
	con_file.write("------------------------<BR>\n")
	if data['time']['tstart']:	
		con_file.write("tstart = " + data['time']['tstart'] + "<BR>\n")
	if data['time']['tend']:	
		con_file.write("tend = " + data['time']['tend'] + "<BR>\n")
	if data['time']['tend'] and data['time']['tstart']:
		con_file.write("# of time steps = " + str(timestp) + "<BR>\n")
	con_file.write('<BR\n>')
	con_file.write("Parameters<BR>\n")
	con_file.write("------------------------<BR>\n")
	if data['parameters']['ice_limit']:
		con_file.write("ice_limit = " + data['parameters']['ice_limit'] + "<BR>\n")
	if data['parameters']['flow_factor']:
		con_file.write("flow_factor = " + data['parameters']['flow_factor'] + "<BR>\n")
	con_file.write('<BR\n>')
	con_file.write("Options<BR>\n")
	con_file.write("------------------------<BR>\n")
	if data['options']['flow_law']:
		con_file.write("flow_law = " + data['options']['flow_law'] + "<BR>\n")
	if data['options']['evolution']:
		con_file.write("evolution = " + data['options']['evolution'] + "<BR>\n")
	if data['options']['temperature']:
		con_file.write("temperature = " + data['options']['temperature'] + "<BR>\n")
	con_file.write('<BR\n>')
	con_file.write('</HTML>\n')
	con_file.close()

# grabbing data from the job output file from the production run

if options.gis_prod:

	out_file = open(target_html + '/GIS-out-diag.html', 'w')

	out_file.write('<HTML>\n')
        out_file.write('<TITLE>Production Job Output Diagnostics</TITLE>\n') 

	procttl, nonlist, avg2, out_flag, nd_name, ld_name = VV_outprocess.jobprocess(gis_output,'gis5km')

#	if error_flag == 1:
#		out_file.write('<FONT COLOR="purple"><H1>Model run incomplete, pick a new job output file for diagnostics!</H1></FONT>')
	if out_flag == 1:
		out_file.write('<FONT COLOR="red"><H2>Job Output Diagnostics</H2></FONT>')
	else:
		out_file.write('<H2>Job Output Diagnostics</H2>')
#	print procttl

# create iteration plots for proudction simulation

	data_script=options.ncl_path + "/solver_gis.ncl"

	plot_gis_data = "ncl 'nfile=\"" + options.data_path + "" + nd_name + "\"'" + ' ' + \
                     "'lfile=\"" + options.data_path + "" + ld_name + "\"'" + ' ' + \
                     "'PNG=\"" + options.ncl_path + "/gis5km_iter\"'" + ' ' + \
                    data_script + ' ' + "1> /dev/null"
	print options.data_path

	try:
		output = subprocess.call(plot_gis_data, shell=True)
	except:
		print "error formatting iteration plot of production run"
		raise

#transferring thickness plot to www location

	if (options.ncl_path + '/gis5km_iter.png'):
	        iterpic = "mv -f " + options.ncl_path + "/gis5km_iter.png" + " " + target_html + "/"
	        try:
	                output = subprocess.call(iterpic, shell=True)
	        except:
	                print "error moving iter png file"
	                raise

		out_file.write('<TABLE>\n')
		out_file.write('<TR>\n')
		out_file.write('<H4>Iteration Count for Nonlinear and Linear Solver</H4>\n')
		out_file.write('<OBJECT data="gis5km_iter.png" type="image/png" width="1300" height="800" hspace=10 align=left alt="Solver Plots">\n')
		out_file.write('</OBJECT>\n')
		out_file.write('<TR>\n')
		out_file.write('<BR>\n')
		out_file.write('</TABLE>\n')

		out_file.write('<BR>\n')
		out_file.write("Number of Processors = " + str(procttl[-1]) + "<BR>\n")
		out_file.write("Number of Nonlinear Iterations = ")
		for item in nonlist:
			out_file.write(str(item) + ", ")
		out_file.write('<BR>\n')
		if out_flag == 1:
			out_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
		out_file.write("Average Number of Linear Iterations per Time-Step = ")
		for item in avg2:
			out_file.write(str(item) + ", ")
		out_file.write('<BR>\n')
		out_file.write('</HTML>\n')
		out_file.close()

# plot production run output for comparison to the benchmark

	if stock_nc:
		plot_file = open(target_html + '/GIS-plot-diag.html', 'w')

		plot_file.write('<HTML>\n')
       		plot_file.write('<TITLE>Thickness</TITLE>\n')
		plot_file.write('<H2>Thickness Plot</H2>')
		plot_file.write('<TABLE>\n')
		plot_file.write('<TR>\n')
		plot_file.write('<H4>a) Benchmark Ice Thickness</H4>\n')
		plot_file.write('<H4>b) Simulation Ice Thickness</H4>\n')
		plot_file.write('<H4>c) Difference from Benchmark </H4>\n')
		plot_file.write('<OBJECT data="gis5km_thk.png" type="image/png" width="1300" height="800" hspace=10 align=left alt="Thickness Plots">\n')
		plot_file.write('</OBJECT>\n')
		plot_file.write('<TR>\n')
		plot_file.write('<BR>\n')
		plot_file.write('</TABLE>\n')
		plot_file.write('</HTML>\n')
		plot_file.close()

file.write('</BODY>\n')
file.write('</HTML>\n')
file.close()

print "LIVV Completed. Go to " + target_html + "/GIS-main-diag.html to view results"
