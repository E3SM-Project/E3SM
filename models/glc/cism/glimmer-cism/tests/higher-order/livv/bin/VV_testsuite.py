#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_outprocess
import VV_dome30details
import VV_shelfdetails
import VV_ismip
import VV_gis10details

def web(descript_file,test_file, \
	dome30d_file,dome30d_case,dome30d_plot,dome30e_file,dome30e_case,dome30e_plot, \
	circ_file,circ_case,circ_plot,conf_file,conf_case,conf_plot, \
	ishoma80_file,ishoma80_case,ishoma80_plot,ishomc80_file,ishomc80_case,ishomc80_plot, \
	gis10_file,gis10_case,gis10_plot,job_path,ncl_path,data_path,html_path):  

# using data, fill the web page with info about the cases
	test_file.write('<HTML>\n')
	test_file.write('<TITLE>Test Suite Diagnostics</TITLE>\n')
	test_file.write('<FONT COLOR="green"><H2>Test Suite Diagnostics</H2></FONT>')

# link to descript_file about the test cases
	test_file.write('<TH ALIGN=LEFT><A HREF="test_descript.html">Test Suite Descriptions</A>\n')
	test_file.write('<BR>\n')

# Diagnostic Dome 30 stats
	test_file.write('<H3>Diagnostic Dome 30 Test:</H3>')

# put something here to flag BFB results and no need to do any more calculations
	flag_to_plot_dome30d = 1
	if flag_to_plot_dome30d:

# link to dome30d_file with descriptions about the test cases
		test_file.write('<TH ALIGN=LEFT><A HREF="dome30d_details.html">Diagnostic Dome 30 Velocity Solver Details</A>\n')
		test_file.write('<BR>\n')
        	VV_dome30details.ddetails(dome30d_file,job_path,ncl_path,data_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="dome30d_case.html">Diagnostic Dome 30 Case Details</A>\n')
		test_file.write('<BR>\n')
                confpath = job_path + '/dome30/diagnostic/dome.30.JFNK.trilinos.config.1'
                benchpath = job_path + '/dome30/diagnostic/bench/dome.30.JFNK.trilinos.config.1'
        	VV_dome30details.case(dome30d_case,confpath,benchpath,ncl_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="dome30d_plot.html">Diagnostic Dome 30 Plots</A>\n')
		test_file.write('<BR>\n')
        	VV_dome30details.dplot(dome30d_plot,job_path,ncl_path,html_path)

# Evolving Dome 30 stats
	test_file.write('<H3>Evolving Dome 30 Test:</H3>')

# put something here to flag BFB results and no need to do any more calculations
	flag_to_plot_dome30e = 1
	if flag_to_plot_dome30e:

# link to dome30e_file with descriptions about the test cases
		test_file.write('<TH ALIGN=LEFT><A HREF="dome30e_details.html">Evolving Dome 30 Velocity Solver Details</A>\n')
		test_file.write('<BR>\n')
        	VV_dome30details.edetails(dome30e_file,job_path,ncl_path,data_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="dome30e_case.html">Evolving Dome 30 Case Details</A>\n')
		test_file.write('<BR>\n')
                confpath = job_path + '/dome30/evolving/dome.30.JFNK.trilinos.config.15'
                benchpath = job_path + '/dome30/evolving/bench/dome.30.JFNK.trilinos.config.15'
        	VV_dome30details.case(dome30e_case,confpath,benchpath,ncl_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="dome30e_plot.html">Evolving Dome 30 Plots</A>\n')
		test_file.write('<BR>\n')
        	VV_dome30details.eplot(dome30e_plot,job_path,ncl_path,html_path)

# Circular Shelf stats
	test_file.write('<H3>Circular Shelf Test</H3>')

# put something here to flag BFB results and no need to do any more calculations
	flag_to_plot_circ = 1
	if flag_to_plot_circ:

		test_file.write('<TH ALIGN=LEFT><A HREF="circ_details.html">Circular Shelf Velocity Solver Details</A>\n')
		test_file.write('<BR>\n')
        	VV_shelfdetails.circdetails(circ_file,job_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="circ_case.html">Circular Shelf Case Details</A>\n')
		test_file.write('<BR>\n')
                confpath = job_path + '/circular-shelf/circular-shelf.JFNK.config'
                benchpath = job_path + '/circular-shelf/bench/circular-shelf.JFNK.config'
        	VV_dome30details.case(circ_case,confpath,benchpath,ncl_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="circ_plot.html">Circular Shelf Plots</A>\n')
		test_file.write('<BR>\n')
        	VV_shelfdetails.circplot(circ_plot,job_path,ncl_path,html_path)

# Confined Shelf stats
	test_file.write('<H3>Confined Shelf Test</H3>')

# put something here to flag BFB results and no need to do any more calculations
	flag_to_plot_conf = 1
	if flag_to_plot_conf:

		test_file.write('<TH ALIGN=LEFT><A HREF="conf_details.html">Confined Shelf Velocity Solver Details</A>\n')
		test_file.write('<BR>\n')
        	VV_shelfdetails.confdetails(conf_file,job_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="conf_case.html">Confined Shelf Case Details</A>\n')
		test_file.write('<BR>\n')
                confpath = job_path + '/confined-shelf/confined-shelf.JFNK.config'
                benchpath = job_path + '/confined-shelf/bench/confined-shelf.JFNK.config'
        	VV_dome30details.case(conf_case,confpath,benchpath,ncl_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="conf_plot.html">Confined Shelf Plot</A>\n')
		test_file.write('<BR>\n')
        	VV_shelfdetails.confplot(conf_plot,job_path,ncl_path,html_path)

# ISMIP HOM A stats
	test_file.write('<H3>ISMIP HOM A 80km Test</H3>')

# put something here to flag BFB results and no need to do any more calculations
	flag_to_plot_iha = 1
	if flag_to_plot_iha:

		test_file.write('<TH ALIGN=LEFT><A HREF="ishoma80_details.html">ISMIP HOM A 80km Velocity Solver Details</A>\n')
		test_file.write('<BR>\n')
        	VV_ismip.a80details(ishoma80_file,job_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="ishoma80_case.html">ISMIP HOM A 80km Case Details</A>\n')
		test_file.write('<BR>\n')
                confpath = job_path + '/ismip-hom-a/80km/ishom.a.80km.JFNK.trilinos.config'
                benchpath = job_path + '/ismip-hom-a/80km/bench/ishom.a.80km.JFNK.trilinos.config'
        	VV_dome30details.case(ishoma80_case,confpath,benchpath,ncl_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="ishoma80_plot.html">ISMIP HOM A 80km Plots</A>\n')
		test_file.write('<BR>\n')
        	VV_ismip.a80plot(ishoma80_plot,job_path,ncl_path,html_path)

# ISMIP HOM C stats
	test_file.write('<H3>ISMIP HOM C 80km Test</H3>')

# put something here to flag BFB results and no need to do any more calculations
	flag_to_plot_ihc = 1
	if flag_to_plot_ihc:

		test_file.write('<TH ALIGN=LEFT><A HREF="ishomc80_details.html">ISMIP HOM C 80km Velocity Solver Details</A>\n')
		test_file.write('<BR>\n')
        	VV_ismip.c80details(ishomc80_file,job_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="ishomc80_case.html">ISMIP HOM C 80km Case Details</A>\n')
		test_file.write('<BR>\n')
                confpath = job_path + '/ismip-hom-c/80km/ishom.c.80km.JFNK.trilinos.config'
                benchpath = job_path + '/ismip-hom-c/80km/bench/ishom.c.80km.JFNK.trilinos.config'
        	VV_dome30details.case(ishomc80_case,confpath,benchpath,ncl_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="ishomc80_plot.html">ISMIP HOM C 80km Plots</A>\n')
		test_file.write('<BR>\n')
        	VV_ismip.c80plot(ishomc80_plot,job_path,ncl_path,html_path)

# GIS 10km stats

	test_file.write('<H3>GIS 10km Tests</H3>')

# put something here to flag BFB results between data and bench and pgi and gnu etc.
	flag_to_plot_gis10km = 1
	if flag_to_plot_gis10km:

		test_file.write('<TH ALIGN=LEFT><A HREF="gis10_details.html">GIS 10km Velocity Solver Details</A>\n')
		test_file.write('<BR>\n')
        	VV_gis10details.details(gis10_file,job_path,ncl_path,data_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="gis10_case.html">GIS 10km Case Details</A>\n')
		test_file.write('<BR>\n')
                confpath = job_path + '/gis_10km/gis_10km.JFNK.trilinos.10.config'
                benchpath = job_path + '/gis_10km/bench/gis_10km.JFNK.trilinos.10.config'
        	VV_dome30details.case(gis10_case,confpath,benchpath,ncl_path,html_path)

		test_file.write('<TH ALIGN=LEFT><A HREF="gis10_plot.html">GIS 10km Plots</A>\n')
		test_file.write('<BR>\n')
        	VV_gis10details.gis10_plot(gis10_plot,job_path,ncl_path,html_path)

	test_file.write('</HTML>\n')
	test_file.close()

#	descript_file = open(options.html_path + '/test_descript.html', 'w')
	descript_file.write('<HTML>\n')
	descript_file.write('<TITLE>Descriptions about the Test Suite</TITLE>\n')
	descript_file.write('<H2>Test Suite Details</H2>')
	descript_file.write('<BR>\n')
	descript_file.write('The Diagnostic Dome 30 test case \n')
	descript_file.write('<BR>\n')
	descript_file.write('  Attributes \n')
	descript_file.write('<BR>\n')
	descript_file.write('  What does it test? \n')
	descript_file.write('<BR><BR>\n')
	descript_file.write('The Evolving Dome 30 test case \n')
	descript_file.write('<BR>\n')
	descript_file.write('  Attributes \n')
	descript_file.write('<BR>\n')
	descript_file.write('  What does it test? \n')
	descript_file.write('<BR><BR>\n')
	descript_file.write('The Circular Shelf test case \n')
	descript_file.write('<BR>\n')
	descript_file.write('  Attributes \n')
	descript_file.write('<BR>\n')
	descript_file.write('  What does it test? \n')
	descript_file.write('<BR><BR>\n')
	descript_file.write('The Confined Shelf test case \n')
	descript_file.write('<BR>\n')
	descript_file.write('  Attributes \n')
	descript_file.write('<BR>\n')
	descript_file.write('  What does it test? \n')
	descript_file.write('<BR><BR>\n')
	descript_file.write('The ISMIP HOM A 80km test case \n')
	descript_file.write('<BR>\n')
	descript_file.write('  Attributes \n')
	descript_file.write('<BR>\n')
	descript_file.write('  What does it test? \n')
	descript_file.write('<BR><BR>\n')
	descript_file.write('The ISMIP HOM C 80km test case \n')
	descript_file.write('<BR>\n')
	descript_file.write('  Attributes \n')
	descript_file.write('<BR>\n')
	descript_file.write('  What does it test? \n')
	descript_file.write('<BR><BR>\n')
	descript_file.write('The Greenland Ice Sheet 10km test case \n')
	descript_file.write('<BR>\n')
	descript_file.write('  Attributes \n')
	descript_file.write('<BR>\n')
	descript_file.write('  What does it test? \n')
	descript_file.write('<BR><BR>\n')
	descript_file.write('</HTML>\n')
	descript_file.close()

#	dome30_file.write('</HTML>\n')
#	dome30_file.close()
