#!/usr/bin/env python3

import netcdf4_functions as nffun
import os, sys, csv, time, math, numpy
from optparse import OptionParser

#Create, run and process a CLM/ELM model ensemble member
#  given specified case and parameters (see parm_list and parm_data files)
#  Parent case must be pre-built and all namelists in run directory.
#  Post-processing calcuates normalized sum of squared errors (SSE) given
#  data constraints specified in "constraints" directory"
#  DMRicciuto 12/1/2015
#
#  Note:  This will only work for single-point CLM/ELM compiled with MPI_SERIAL

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
parser.add_option("--finidat_thiscase", dest='fincase', default=False, \
                  action="store_true", help = 'Use this case for finidat')
parser.add_option('--spin_cycle', dest='spin_cycle', default=0, \
                  help = 'Number of years in spinup cycle')
parser.add_option('--1850_landuse', dest='nolanduse', default=False, \
                  action='store_true', help = '1850 land use (no dynamics)')
parser.add_option('--1850_co2', dest='noco2', default=False, \
                  action='store_true', help = 'use 1850 CO2')
parser.add_option('--1850_ndep', dest='nondep', default=False, \
                  action='store_true', help = 'use 1850 NDep')
parser.add_option('--suffix', dest='suffix', default='', \
                  help = 'use 1850 NDep')
parser.add_option('--machine', dest='machine', default='cades', \
                  help = 'machine')
parser.add_option('--warming', dest='warming', default='0.0', \
                  help = 'warming level to apply')
(options, args) = parser.parse_args()


casename = options.casename

#create directory from original case 
orig_dir = str(os.path.abspath(options.runroot)+'/'+casename+'/run')
new_dir  = orig_dir.replace(options.site_orig, options.site_new)
if (options.suffix != ''):
  new_dir = new_dir.replace(casename,casename+'_'+options.suffix)
print 'Copying from '+orig_dir+' to \n'+new_dir
if (new_dir == orig_dir):
  print 'Error:  New and old directories are the same.  Exiting'
  sys.exit(1)

#copy files to new directory
os.system('mkdir -p '+new_dir+'/timing/checkpoints')
os.system('rm  '+new_dir+'/*.nc')
os.system('rm -f '+new_dir+'/*.nc.orig')
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
            elif ('drv_in' in f and 'lnd_ntasks' in s):
                np = int(s.split()[2])
                s_out = s
            elif ('drv_in' in f and 'restart_n' in s and int(options.nyears) > 0):
                s_out = ' restart_n = '+str(options.nyears)+'\n'
            elif ('lnd_in' in f and 'finidat =' in s and int(options.finyr) > 0):
                year_orig = str((s.split('.')[-2:-1])[0])[0:4]
                year_new = str(10000+int(options.finyr))[1:]
                s_out = s.replace('.clm2.r.'+year_orig, '.clm2.r.'+year_new)
                s_out = s_out.replace(options.site_orig, options.site_new)
                if (options.suffix != ''):
                  s_out = s_out.replace(casename,casename+'_'+options.suffix)
            elif ('lnd_in' in f and "hist_nhtfrq =" in s and \
                        int(options.spin_cycle) > 0):
                nhtfrq = str(-8760*int(options.spin_cycle))
                if ('ad_spinup' in new_dir):
                    s_out = 'hist_nhtfrq = '+nhtfrq+','+nhtfrq+'\n'
                else:
                    s_out = 'hist_nhtfrq = '+nhtfrq+'\n'
            elif ('lnd_in' in f and "flanduse_timeseries =" in s and \
                     options.nolanduse == True):
               s_out = " flanduse_timeseries = ''\n"
            elif ('lnd_in' in f and "do_transient_pfts" in s and \
                     options.nolanduse == True):
               s_out = "do_transient_pfts = .false.\n"
            elif ('lnd_in' in f and "do_harvest" in s and \
                     options.nolanduse == True):
               s_out = "do_harvest = .false.\n"
            elif ('lnd_in' in f and "co2_file =" in s and \
                     options.noco2 == True):
               s_out = s.replace('.nc','_CON.nc')
            elif ('lnd_in' in f and "stream_fldfilename_ndep" in s and \
                     options.nondep == True):
               s_out = s.replace('.nc','_CON.nc')
            elif ('diri =' not in s):
                s_out = s.replace(options.site_orig, options.site_new)
                if (options.suffix != ''):
                  s_out = s_out.replace(casename,casename+'_'+options.suffix)
            elif ('diri' in s and 'lnd' in f):
                  exedir = s.split()[2][1:-4]
                  print exedir
                  s_out = s
            else:
                s_out = s
            myoutput.write(s_out)
        myoutput.close()
        myinput.close()
        os.system(' mv '+new_dir+'/'+f+'.tmp '+new_dir+'/'+f)
#Assume makepointdata has been run to generate surface and domain data
if (options.site_orig == options.site_new and os.path.exists(orig_dir+'/surfdata.nc')):
  os.system('cp '+orig_dir+'/surfdata.nc '+new_dir)
else:
  os.system('cp temp/surfdata.nc '+new_dir)
if (options.site_orig == options.site_new and os.path.exists(orig_dir+'/domain.nc')):
  os.system('cp '+orig_dir+'/domain.nc '+new_dir)
else:
  os.system('cp temp/domain.nc '+new_dir)


os.system('pwd')
if ('20TR' in options.casename and not options.nolanduse):
   os.system('cp temp/*pftdyn*.nc '+new_dir)


#if a global file exists, modify
if (os.path.exists('temp/global_'+options.casename+'_0.pbs') and options.suffix != ''):
  file_in = open('temp/global_'+options.casename+'_0.pbs','r')
  file_out =open('temp/global_'+options.casename+'_'+options.suffix+'.pbs','w')
  mpicmd = 'mpirun -np '+str(np)
  if ('cades' in options.machine):
    mpicmd = '/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/1.10.3/centos7.2_gnu5.3.0/bin/mpirun'+ \
            ' -np '+str(np)+' --hostfile $PBS_NODEFILE ' 
  elif (('titan' in options.machine or 'eos' in options.machine) and int(options.ninst) == 1):
    mpicmd = 'aprun -n '+str(np)
  elif ('cori' in options.machine or 'edison' in options.machine):
    mpicmd = 'srun -n '+str(np)

  for s in file_in:
    if "#" in s:
      file_out.write(s)
  file_out.write('\n\n')
  file_out.write('cd '+new_dir+'\n')
  file_out.write(mpicmd+' '+exedir+'/e3sm.exe\n')
  file_in.close()
  file_out.close()
  print "Submitting the job:"
  if ('cori' in options.machine or 'edison' in options.machine):
    os.system('sbatch temp/global_'+options.casename+'_'+options.suffix+'.pbs')
  else:
    os.system('qsub temp/global_'+options.casename+'_'+options.suffix+'.pbs')

#if an ensemble script exists, make a copy for this site
if (options.site_orig != ''):
  case_prefix = options.casename.split(options.site_orig)[0][:-1]
  ens_fname = 'scripts/'+case_prefix+'/ensemble_run_'+options.casename.replace(options.site_orig, \
                 options.site_new)+'.pbs'

  if os.path.exists('scripts/'+case_prefix+'/ensemble_run_'+options.casename+'.pbs'):
    os.system('cp scripts/'+case_prefix+'/ensemble_run_'+options.casename+'.pbs '+ens_fname)
    myfile = open(ens_fname,'r')
    myfile_out = open(ens_fname+'_temp','w')

    for s in myfile:
      args = s.split(' ')
      if ('--case' in args):
        caseind = args.index('--case')
        print caseind
        args[caseind+1] = args[caseind+1].replace(options.site_orig,options.site_new)
        siteind = args.index('--site')
        args[siteind+1] = args[siteind+1].replace(options.site_orig,options.site_new)
        sout = " ".join(args)
        myfile_out.write(sout)
      else:
        myfile_out.write(s.replace(options.site_orig,options.site_new))
    myfile.close()
    myfile_out.close()
    os.system('mv '+ens_fname+'_temp '+ens_fname)  

