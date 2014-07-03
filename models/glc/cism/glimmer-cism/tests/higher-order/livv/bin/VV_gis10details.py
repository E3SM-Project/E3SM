#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess
import collections
import VV_outprocess
import VV_utilities
import VV_checks

def details(solver_file,job_path,ncl_path,data_path,target_html):  # using data, fill the web page with info
        
        failedt_list = []

	solver_file.write('<HTML>\n')
	solver_file.write('<H3>GIS 10km Iteration Count Details:</H3>')
	solver_file.write('<H4>Eventually published in plot form</H4>')
	solver_file.write('<BR> \n')

# JFNK multiple procs
        
# Failure checking
        failedt1 = VV_checks.failcheck(job_path, '/gis_10km/data/gis_10km.JFNK.trilinos.10.gnu.out')
        failedt_list.append(failedt1)

        procttl_gis10d, nonlist_gis10d, avg2_gis10d, out_flag_gis10d, ndg102_name, ldg102_name = VV_outprocess.jobprocess(job_path + '/gis_10km/data/gis_10km.JFNK.trilinos.10.gnu.out','gis10km2')

        procttl_gis10b, nonlist_gis10b, avg2_gis10b, out_flag_gis10b, ndg102b_name, ldg102b_name = VV_outprocess.jobprocess(job_path + '/bench/gis_10km/data/gis_10km.JFNK.trilinos.10.gnu.out','gis10km2b')


# create iteration plots for proudction simulation

#        data_script=ncl_path + "/solver_gis10.ncl" 

#        plot_gis10_data = "ncl 'nfile=\"" + data_path + "" + ndg102_name + "\"'" + ' ' + \
#                     "'lfile=\"" + data_path + "" + ldg102_name + "\"'" + ' ' + \
#                     "'nbfile=\"" + data_path + "" + ndg102b_name + "\"'" + ' ' + \
#                     "'lbfile=\"" + data_path + "" + ldg102b_name + "\"'" + ' ' + \
#                     "'PNG=\"" + ncl_path + "/gis10km_iter\"'" + ' ' + \
#                    data_script + ' ' + "1> /dev/null"

#        try:
#                output = subprocess.call(plot_gis10_data, shell=True)
#        except:
#                print "error formatting iteration plot of gis10km run"
#                raise

#transferring iteration plot to www location

#        if (ncl_path + '/gis10km_iter.png'):
#                iterpic = "mv -f " + ncl_path + "/gis10km_iter.png" + " " + target_html + "/"
#                try:
#                        output = subprocess.call(iterpic, shell=True)
#                except:
#                        print "error moving iter png file"
#                        raise

#        solver_file.write('<TABLE>\n')
#        solver_file.write('<TR>\n') 
#        solver_file.write('<H4>Iteration Count for Nonlinear and Linear Solver</H4>\n')
#        solver_file.write('<OBJECT data="gis10km_iter.png" type="image/png" width="1300" height="800" hspace=10 align=left alt="Solver Plots">\n')
#        solver_file.write('</OBJECT>\n')
#        solver_file.write('<TR>\n')
#        solver_file.write('<BR>\n')
#        solver_file.write('</TABLE>\n')

# also present data in list form
#        if fatal_flag == 1:
#                solver_file.write('<H4>Fatal Error Found in Output File</H4>\n')
        solver_file.write('<H4>New Run: gis_10km.JFNK.trilinos.10.gnu.out</H4>')
        solver_file.write("Number of Processors = " + str(procttl_gis10d[-1]) + "<BR>\n")
        solver_file.write("Number of Nonlinear Iterations = ")
        VV_utilities.format(solver_file, nonlist_gis10d)
        solver_file.write('<BR>\n')
        if out_flag_gis10d == 1:
                solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        solver_file.write("Average Number of Linear Iterations per Time-Step = ")
        VV_utilities.format(solver_file, avg2_gis10d)
        solver_file.write('<BR> \n')

        solver_file.write('<H4>Benchmark Run: gis_10km.JFNK.trilinos.10.gnu.out</H4>')
        solver_file.write("Number of Processors = " + str(procttl_gis10b[-1]) + "<BR>\n")
        solver_file.write("Number of Nonlinear Iterations = ")
        VV_utilities.format(solver_file, nonlist_gis10b)
        solver_file.write('<BR>\n')
        if out_flag_gis10b == 1:
                solver_file.write('<FONT COLOR="red">***TIME STEP(S) WHICH FAILED TO CONVERGE</FONT> <BR>\n')
        solver_file.write("Average Number of Linear Iterations per Time-Step = ")
        VV_utilities.format(solver_file, avg2_gis10b)
        solver_file.write('<BR> \n')

        solver_file.write('</HTML>\n')
        solver_file.close()
    
        if 1 in failedt_list:
            failedt = 1
        else:
            failedt = 0

        return failedt

def gis10_plot(plot_file,job_path,ncl_path,html_path):  # using data, fill the web page with info

        tmpath = job_path + '/gis_10km/data/gis_10km.seacism.nc'
        if VV_utilities.emptycheck(tmpath) == 0:
                return

        plot_file.write('<HTML>\n')
        plot_file.write('<H3>GIS 10km Plot Details:</H3>')
        gis10km_plotfile=''+ ncl_path + '/gis10km.ncl'
        stock='STOCK = addfile(\"'+ job_path + '/bench/gis_10km/data/gis_10km.seacism.nc\", \"r\")'
        VAR  ='VAR = addfile(\"' + job_path + '/gis_10km/data/gis_10km.seacism.nc\", \"r\")'
        png  = 'PNG = "' + ncl_path + '/gis10km"'
        plot_gis10km = "ncl '" + VAR + "' '" + stock + "' '" + png + "' " + gis10km_plotfile

        try:
                output = subprocess.call(plot_gis10km, shell=True)
        except:
                print "error creating ncl gis10km plot"
                raise

# transferring conf pic to www file

        if (ncl_path + '/gis10km.png'):
                gispic = "mv -f " + ncl_path + "/gis10km.png" + " " + html_path + "/"
                try:
                        output = subprocess.call(gispic, shell=True)
                except:
                        print "error moving gis10km png file to www directory"
                        sys.exit(1)
                        raise

        plot_file.write('<HTML>\n')
        plot_file.write('<TITLE>GIS 10km Test Case </TITLE>\n')
        plot_file.write('<TABLE>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<H4>Difference from benchmark for a range of processor counts for a range of variables</H4>\n')
        plot_file.write('<OBJECT data="gis10km.png" type="image/png" width="1100" height="800" hspace=10 align=left alt="GIS 10km Plots PNG">\n')
        plot_file.write('</OBJECT>\n')
        plot_file.write('<TR>\n')
        plot_file.write('<BR>\n')
        plot_file.write('</TABLE>\n')

        plot_file.write('</HTML>\n')
        plot_file.close()

