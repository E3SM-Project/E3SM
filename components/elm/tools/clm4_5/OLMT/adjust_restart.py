#!/usr/bin/env/python
import os, sys, csv, time, math
from optparse import OptionParser
import netcdf4_functions as nffun
import glob

parser = OptionParser()

parser.add_option('--rundir', dest='rundir', default="", \
                    help = 'location of run directory')
parser.add_option('--casename', dest='casename', default="", \
                    help='Name of case')
parser.add_option('--restart_year', dest='restart_year', default="", \
                    help='Year of restart file to modify')
parser.add_option('--BGC', dest='bgc', default=False, action="store_true", \
                    help='Flag to set for BGC compsets')
parser.add_option('--harvest', dest='harvest', default=False, action="store_true", \
                    help='Add a 95% above ground harvest')

(options,args)=parser.parse_args()

#The purpose of this code is to replace the restart value with the average value of the last met cycle
#This will bring the accelerated pools more into line with the expected results.

casename = options.casename
if (options.restart_year == ''):
        #if restart_year not provided, take the last existing restart file
	restart_file = glob.glob(options.rundir+'/'+casename+'.clm2.r.*.nc')
        if (len(restart_file) > 1):
    	   restart_file_last = restart_file[-1]
        else:
           restart_file_last = restart_file[0]
	year = int(restart_file_last[-19:-15])
else:
	year = int(options.restart_year)

if ('BGC' in casename):
	options.bgc = True

fname_restart = options.rundir+'/'+casename+'.clm2.r.'+str(10000+year)[1:]+ \
    '-01-01-00000.nc'
fname_hist    = options.rundir+'/'+casename+'.clm2.h1.'+str(10000+year)[1:]+ \
    '-01-01-00000.nc'

#save original restart file
if (os.path.isfile(fname_restart+'.orig')):
  os.system('cp '+fname_restart+'.orig '+fname_restart)
else:
  os.system('cp '+fname_restart+' '+fname_restart+'.orig')
os.system('chmod a-w '+fname_restart+'.orig')

#Accelerated pools:  Coarse woody debris, dead stem, coarse root, soil3 and soil4 C and N
var_names   = ['DEADSTEMC', 'DEADSTEMN', 'DEADSTEMP', 'DEADCROOTC', 'DEADCROOTN', 'DEADCROOTP']
if (options.bgc):
  var_names2d = ['CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL3C_vr', 'SOIL3N_vr', 'SOIL3P_vr', 'SOIL2C_vr', \
                   'SOIL2N_vr', 'SOIL2P_vr']
else:
  var_names2d = ['CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL4C_vr', 'SOIL4N_vr', 'SOIL4P_vr', 'SOIL3C_vr', \
                   'SOIL3N_vr', 'SOIL3P_vr']
if (options.harvest):
  var_names_harvest = ['DEADSTEMC','DEADSTEMN','DEADSTEMP','LIVESTEMN','LIVESTEMC','LIVESTEMP', \
                       'LEAFC', 'LEAFN', 'LEAFP', 'DEADSTEMC_STORAGE', 'DEADSTEMN_STORAGE', \
	               'DEADSTEMP_STORAGE', 'LIVESTEMC_STORAGE', 'LIVESTEMN_STORAGE', \
                       'LIVESTEMP_STORAGE', 'LEAFC_STORAGE', 'LEAFN_STORAGE', 'LEAFP_STORAGE', \
                       'FROOTC_STORAGE', 'FROOTC', 'FROOTN_STORAGE', 'FROOTN', 'FROOTP_STORAGE', \
                       'FROOTP', 'LIVECROOTC', 'LIVECROOTC_STORAGE', 'LIVECROOTN', 'LIVECROOTN_STORAGE', \
                       'LIVECROOTP', 'LIVECROOTP_STORAGE', 'DEADCROOTC', 'DEADCROOTC_STORAGE', \
                       'DEADCROOTN', 'DEADCROOTN_STORAGE', 'DEADCROOTP', 'DEADCROOTP_STORAGE', 'XSMRPOOL', \
                       'XSMRPOOL_RECOVER']

if (options.harvest):
  for v in range(0,len(var_names_harvest)):	
    rest_vals = nffun.getvar(fname_restart, var_names_harvest[v].lower())
    #Loop through all valid values in the restart file.
    n_rest = len(rest_vals)
    for i in range(0,n_rest):
      rest_vals[i] = rest_vals[i]*0.05
      if (var_names_harvest[v] == 'LEAFC'):
        rest_vals[i] = 0.33/0.03
        print rest_vals[i]
      elif (var_names_harvest[v] == 'FROOTC'):
        rest_vals[i] = 0.33/0.03*0.666
      elif (var_names_harvest[v] == 'LEAFN'):
        rest_vals[i] = 0.33/0.03/25.0
      elif (var_names_harvest[v] == 'FROOTN'):
        rest_vals[i] = 0.33/0.03/42.0
    ierr = nffun.putvar(fname_restart, var_names_harvest[v].lower(), rest_vals)
else:
  for v in range(0,len(var_names)):
    hist_vals = nffun.getvar(fname_hist, var_names[v])
    rest_vals = nffun.getvar(fname_restart, var_names[v].lower())
    #Loop through all valid values in the restart file.
    n_rest = len(rest_vals)
    for i in range(0,n_rest):
      if (float(rest_vals[i]) > 0.0 and float(hist_vals[0][i]) < 1.0e10 and float(hist_vals[0][i]) > 0.001):
        rest_vals[i] = hist_vals[0][i]
    ierr = nffun.putvar(fname_restart, var_names[v].lower(), rest_vals)

  #get a single, non-depth dependent variable to count # of columns
  rest_vals = nffun.getvar(fname_restart, 'fpg')
  n_rest = len(rest_vals) 
  for v in range(0,len(var_names2d)):
    hist_vals = nffun.getvar(fname_hist, var_names2d[v])
    rest_vals = nffun.getvar(fname_restart, var_names2d[v].lower())
    #Loop through all valid values in the restart file.
    for i in range(0,n_rest):
      for j in range(0,10):
        if (float(rest_vals[i][j]) > 0.0 and float(hist_vals[0][j][i]) < 1.0e10 and float(hist_vals[0][j][i]) > 1e-10):
          rest_vals[i][j] = hist_vals[0][j][i]
    ierr = nffun.putvar(fname_restart, var_names2d[v].lower(), rest_vals)

#Remove negative Ppool values
#os.system("ncap2 -O -s 'ppool=0.2*npool' "+fname_restart+" "+fname_restart)
#os.system("ncap2 -O -s 'defdim(\"cohort\",5663)' -s 'defdim(\"levurb\",5)' -s " \
#   +"'defdim(\"string_length\",64)' -s 'defdim(\"levtrc\",10)' " + \
#   fname_restart+" "+fname_restart)

