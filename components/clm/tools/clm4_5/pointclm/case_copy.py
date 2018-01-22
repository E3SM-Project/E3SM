#!/usr/bin/env python

import netcdf_functions as nffun
import os, sys, csv, time, math, numpy
from optparse import OptionParser

#Create, run and process a CLM/ALM model ensemble member
#  given specified case and parameters (see parm_list and parm_data files)
#  Parent case must be pre-built and all namelists in run directory.
#  Post-processing calcuates normalized sum of squared errors (SSE) given
#  data constraints specified in "constraints" directory"
#  DMRicciuto 12/1/2015
#
#  Note:  This will only work for single-point CLM/ALM compiled with MPI_SERIAL

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--runroot", dest="runroot", default="../../run", \
                  help="Directory where the run would be created")
parser.add_option("--case_copy", dest="casename", default="", \
                  help="Name of case to copy from")
parser.add_option("--site_orig", dest="site_orig", default='', \
                  help = 'Site being run in original case')
parser.add_option("--site_new", dest="site_new", default='', \
                  help = 'Site to run')
parser.add_option("--nyears", dest="nyears", default=0, \
                  help = 'Number of years to run')
parser.add_option("--finidat_year", dest='finyr', default=0, \
                  help = 'Year for initial data file')
parser.add_option('--spin_cycle', dest='spin_cycle', default=0, \
                  help = 'Number of years in spinup cycle')

(options, args) = parser.parse_args()

casename = options.casename

#create directory from original case 
orig_dir = str(os.path.abspath(options.runroot)+'/'+casename+'/run')
new_dir  = orig_dir.replace(options.site_orig, options.site_new)
print 'Copying from '+orig_dir+' to \n'+new_dir

#copy files to new directory
os.system('mkdir -p '+new_dir+'/timing/checkpoints')
os.system('rm  '+new_dir+'/*.nc')
os.system('cp  '+orig_dir+'/*_in* '+new_dir)
os.system('cp  '+orig_dir+'/*nml '+new_dir)
if (not ('ICB' in casename)):
    os.system('cp  '+orig_dir+'/*stream* '+new_dir)
os.system('cp  '+orig_dir+'/*.rc '+new_dir)
os.system('cp  '+orig_dir+'/*para*.nc '+new_dir)

#Change site name in relevant files

for f in os.listdir(new_dir):
    if (os.path.isfile(new_dir+'/'+f) and (f[-2:] == 'in' or f[-3:] == 'nml' or 'streams' in f)):
        myinput=open(new_dir+'/'+f)
        myoutput=open(new_dir+'/'+f+'.tmp','w')
        for s in myinput:
            if ('drv_in' in f and 'stop_n' in s and int(options.nyears) > 0):
                s_out = ' stop_n = '+str(options.nyears)+'\n'
            elif ('drv_in' in f and 'restart_n' in s and int(options.nyears) > 0):
                s_out = ' restart_n = '+str(options.nyears)+'\n'
            elif ('lnd_in' in f and 'finidat =' in s and int(options.finyr) > 0):
                year_orig = str((s.split('.')[-2:-1])[0])[0:4]
                year_new = str(10000+int(options.finyr))[1:]
                s_out = s.replace('.clm2.r.'+year_orig, '.clm2.r.'+year_new)
                s_out = s_out.replace(options.site_orig, options.site_new)
            elif ('lnd_in' in f and "hist_nhtfrq =" in s and \
                        int(options.spin_cycle) > 0):
                nhtfrq = str(-8760*int(options.spin_cycle))
                if ('ad_spinup' in new_dir):
                    s_out = 'hist_nhtfrq = '+nhtfrq+','+nhtfrq+'\n'
                else:
                    s_out = 'hist_nhtfrq = '+nhtfrq+'\n'
            elif ('diri =' not in s):
                s_out = s.replace(options.site_orig, options.site_new)
            else:
                s_out = s
            myoutput.write(s_out)
        myoutput.close()
        myinput.close()
        os.system(' mv '+new_dir+'/'+f+'.tmp '+new_dir+'/'+f)
#Assume makepointdata has been run to generate surface and domain data
os.system('cp temp/surfdata.nc '+new_dir)
os.system('cp temp/domain.nc '+new_dir)
os.system('pwd')
if ('20TR' in options.casename):
   os.system('cp temp/*pftdyn*.nc '+new_dir)



