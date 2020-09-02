#!/usr/bin/env python

import netcdf4_functions as nffun
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
parser.add_option("--ens_num", dest="ens_num", default=1, \
                  help="Ensemble member number")
parser.add_option("--case", dest="casename", default="", \
                  help="Name of case")
parser.add_option("--ens_file", dest="ens_file", default="", \
                  help="Name of samples file")
parser.add_option("--parm_list", dest="parm_list", default='parm_list', \
                  help = 'File containing list of parameters to vary')
parser.add_option("--cnp", dest="cnp", default = False, action="store_true", \
                  help = 'CNP mode - initialize P pools')
parser.add_option("--site", dest="site", default='parm_list', \
                  help = 'Site name')
(options, args) = parser.parse_args()


parm_names=[]
parm_indices=[]
parm_values=[]
myinput = open(options.parm_list, 'r')
casename = options.casename

#get parameter names and PFT information
pnum=0
for s in myinput:
   pdata = s.split()
   parm_names.append(pdata[0])
   if (pdata[0] == 'co2'):
     pnum_co2 = pnum
   if (len(pdata) == 3):
     parm_indices.append(-1)
   else:
     parm_indices.append(int(pdata[1]))
   pnum=pnum+1
myinput.close()

#get parameter values
if (options.ens_file == ''):
   myinput = open('./parm_data', 'r')
   for s in myinput:    
      parm_values.append(float(s))
   myinput.close()
   os.system('rm ./parm_data')
else:
   myinput = open(options.ens_file, 'r')
   linenum = 1
   for s in myinput:
      if (int(options.ens_num) == linenum):
         parm_values_str = s.split()
         for v in parm_values_str:
            parm_values.append(float(v))
      linenum = linenum+1
   myinput.close()

n_parameters = len(parm_names)
gst=str(100000+int(options.ens_num))

#create ensemble directory from original case 
est = str(100000+int(options.ens_num))
orig_dir = str(os.path.abspath(options.runroot)+'/'+casename+'/run')
ens_dir  = os.path.abspath(options.runroot)+'/UQ/'+casename+'/g'+gst[1:]
		
os.system('mkdir -p '+options.runroot+'/UQ/'+casename+'/g'+gst[1:]+'/timing/checkpoints')
os.system('cp  '+orig_dir+'/*_in* '+ens_dir)
os.system('cp  '+orig_dir+'/*nml '+ens_dir)
if (not ('CB' in casename)):
    os.system('cp  '+orig_dir+'/*stream* '+ens_dir)
os.system('cp  '+orig_dir+'/*.rc '+ens_dir)
os.system('cp  '+orig_dir+'/surf*.nc '+ens_dir)
os.system('cp  '+orig_dir+'/domain*.nc '+ens_dir)
os.system('cp  '+orig_dir+'/*para*.nc '+ens_dir)

#loop through all filenames, change directories in namelists, change parameter values
for f in os.listdir(ens_dir):
    if (os.path.isfile(ens_dir+'/'+f) and (f[-2:] == 'in' or f[-3:] == 'nml' or 'streams' in f)):
        myinput=open(ens_dir+'/'+f)
        myoutput=open(ens_dir+'/'+f+'.tmp','w')
        for s in myinput:
            if ('fates_paramfile' in s):
                paramfile_orig = ((s.split()[2]).strip("'"))
                if (paramfile_orig[0:2] == './'):
                  paramfile_orig = orig_dir+'/'+paramfile_orig[2:]
                paramfile_new  = './fates_params_'+est[1:]+'.nc'
                os.system('cp '+paramfile_orig+' '+ens_dir+'/'+paramfile_new)
                os.system('nccopy -3 '+ens_dir+'/'+paramfile_new+' '+ens_dir+'/'+paramfile_new+'_tmp')
                os.system('mv '+ens_dir+'/'+paramfile_new+'_tmp '+ens_dir+'/'+paramfile_new)
                myoutput.write(" fates_paramfile = '"+paramfile_new+"'\n")
                fates_paramfile = ens_dir+'/fates_params_'+est[1:]+'.nc'
            elif ('paramfile' in s):
                paramfile_orig = ((s.split()[2]).strip("'"))
                if (paramfile_orig[0:2] == './'):
                   paramfile_orig = orig_dir+'/'+paramfile_orig[2:]
                paramfile_new  = './clm_params_'+est[1:]+'.nc'
                os.system('cp '+paramfile_orig+' '+ens_dir+'/'+paramfile_new)
                os.system('nccopy -3 '+ens_dir+'/'+paramfile_new+' '+ens_dir+'/'+paramfile_new+'_tmp')
                os.system('mv '+ens_dir+'/'+paramfile_new+'_tmp '+ens_dir+'/'+paramfile_new)
                myoutput.write(" paramfile = '"+paramfile_new+"'\n")
                pftfile = ens_dir+'/clm_params_'+est[1:]+'.nc'
            elif ('ppmv' in s and 'co2' in parm_names):
                myoutput.write(" co2_ppmv = "+str(parm_values[pnum_co2])+'\n')
            elif ('fsoilordercon' in s):
                CNPfile_orig = ((s.split()[2]).strip("'"))
                if (CNPfile_orig[0:2] == './'):
                   CNPfile_orig  = orig_dir+'/'+CNPfile_orig[2:]
                CNPfile_new  = './CNP_parameters_'+est[1:]+'.nc'
                os.system('cp '+CNPfile_orig+' '+ens_dir+'/'+CNPfile_new)
                os.system('nccopy -3 '+ens_dir+'/'+CNPfile_new+' '+ens_dir+'/'+CNPfile_new+'_tmp')
                os.system('mv '+ens_dir+'/'+CNPfile_new+'_tmp '+ens_dir+'/'+CNPfile_new)
                myoutput.write(" fsoilordercon = '"+CNPfile_new+"'\n")
                CNPfile = ens_dir+'/CNP_parameters_'+est[1:]+'.nc'
            elif ('fsurdat =' in s):
                surffile_orig = ((s.split()[2]).strip("'"))
                if (surffile_orig[0:2] == './'):
                  surffile_orig = orig_dir+'/'+surffile_orig[2:]
                surffile_new = './surfdata_'+est[1:]+'.nc'
                os.system('cp '+surffile_orig+' '+ens_dir+'/'+surffile_new)
                os.system('nccopy -3 '+ens_dir+'/'+surffile_new+' '+ens_dir+'/'+surffile_new+'_tmp')
                os.system('mv '+ens_dir+'/'+surffile_new+'_tmp '+ens_dir+'/'+surffile_new)
                myoutput.write(" fsurdat = '"+surffile_new+"'\n")
                surffile = ens_dir+'/surfdata_'+est[1:]+'.nc'
            elif ('finidat = ' in s):
                finidat_file_orig = ((s.split()[2]).strip("'"))
                if (finidat_file_orig.strip() != ''):
                   finidat_file_new  = ens_dir+'/'+(finidat_file_orig.split('/')[-1:])[0]
                   if (finidat_file_orig[0:2] == './'):
                      finidat_file_orig = orig_dir+'/'+finidat_file_orig[2:]
                   #get finidat files from previous ensemble cases if available
                   if ('1850' in casename and not ('ad_spinup' in casename)): 
                      finidat_file_path = os.path.abspath(options.runroot)+'/UQ/'+casename.replace('1850CNP','1850CN')+'_ad_spinup/g'+gst[1:]
                      if (os.path.exists(finidat_file_path)):
	                  finidat_file_orig = finidat_file_path+'/*.clm2.r.*.nc'
                          os.system('python adjust_restart.py --rundir '+finidat_file_path+' --casename '+ \
                                      casename.replace('1850CNP','1850CN')+'_ad_spinup')
                   if ('20TR' in casename):
                      finidat_file_path = os.path.abspath(options.runroot)+'/UQ/'+casename.replace('20TR','1850')+ \
                                       '/g'+gst[1:]
                      if (os.path.exists(finidat_file_path)):
                          finidat_file_orig = finidat_file_path+'/*.clm2.r.*.nc'
                          os.system('rm '+finidat_file_path+'/*ad_spinup*.clm2.r.*.nc')
                   os.system('cp '+finidat_file_orig+' '+finidat_file_new)
                   myoutput.write(" finidat = '"+finidat_file_new+"'\n")
                else:
                   myoutput.write(s)
            elif ('logfile =' in s):
                os.system('date +%y%m%d-%H%M%S > mytime'+str(options.ens_num))
                mytinput=open('./mytime'+str(options.ens_num),'r')
                for st in mytinput:
                    timestr = st.strip()
                mytinput.close()
                os.system('rm mytime'+str(options.ens_num))
                myoutput.write(s.replace('`date +%y%m%d-%H%M%S`',timestr))
            else:
                myoutput.write(s.replace(orig_dir,ens_dir))
        myoutput.close()
        myinput.close()
        os.system(' mv '+ens_dir+'/'+f+'.tmp '+ens_dir+'/'+f)

pnum = 0
CNP_parms = ['ks_sorption', 'r_desorp', 'r_weather', 'r_adsorp', 'k_s1_biochem', 'smax', 'k_s3_biochem', \
             'r_occlude', 'k_s4_biochem', 'k_s2_biochem']

for p in parm_names:
   if ('INI' in p):
      if ('BGC' in casename):
         scalevars = ['soil3c_vr','soil3n_vr','soil3p_vr']
      else:
         scalevars = ['soil4c_vr','soil4n_vr','soil4p_vr']
      sumvars = ['totsomc','totsomp','totcolc','totcoln','totcolp']
      for v in scalevars:
         myvar = nffun.getvar(finidat_file_new, v)
         myvar = parm_values[pnum] * myvar
         ierr = nffun.putvar(finidat_file_new, v, myvar)
   elif (p == 'lai'):
     myfile = surffile
     param = nffun.getvar(myfile, 'MONTHLY_LAI')
     param[:,:,:,:] = parm_values[pnum]
     ierr = nffun.putvar(myfile, 'MONTHLY_LAI', param)
   elif (p != 'co2'):
      if (p in CNP_parms):
         myfile= CNPfile
      elif ('fates' in p):
         myfile = fates_paramfile
      else:
         myfile = pftfile
      param = nffun.getvar(myfile, p)
      if ('fates_prt_nitr_stoich_p1' in p):
        #this is a 2D parameter.
         param[parm_indices[pnum] % 6:parm_indices[pnum] % 6+2,parm_indices[pnum]/6] = parm_values[pnum]
         param[parm_indices[pnum] % 6:parm_indices[pnum] % 6+4,parm_indices[pnum]/6] = parm_values[pnum]
      elif (parm_indices[pnum] > 0):
         param[parm_indices[pnum]] = parm_values[pnum]
      elif (parm_indices[pnum] == 0):
         param = parm_values[pnum]
      else:
         param[...] = parm_values[pnum]
      ierr = nffun.putvar(myfile, p, param)
      if ('fr_flig' in p):
         param=nffun.getvar(myfile, 'fr_fcel')
         param[...]=1.0-parm_values[pnum]-parm_values[pnum-1]
         ierr = nffun.putvar(myfile, 'fr_fcel', param)
   pnum = pnum+1
