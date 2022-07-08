#!/usr/bin/env python3

import os, sys, csv, time, math, numpy, getpass
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

parser.add_option("--runroot", dest="runroot", default="", \
                  help="Directory where the run would be created")
parser.add_option("--ens_num", dest="ensnum", default=1, \
                  help="Ensemble member number")
parser.add_option("--parm_list", dest="parm_list", default="", \
                  help="File containing parameter names/pfts to modify")
parser.add_option("--parm_data", dest="parm_data", default="", \
                  help="File containing parameter values and ranges")
parser.add_option("--constraints", dest="constraints", default="", \
                  help="Directory containing constraining variables")
parser.add_option("--norun", dest="norun", default=False, action="store_true", \
                  help="Don't run model (use for testing purposes)")
parser.add_option("--machine", dest="machine", default="cades", \
                  help="My machine")
parser.add_option("--casename", dest="casename", default="", \
                  help = "Name of case to run")

(options, args) = parser.parse_args()

#================= netcdf manipulation functions ===============================#
def getvar(fname, varname):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"r",mmap=False)
        var = nffile.variables[varname]
        varvals = var[:].copy()    #works for vector only?
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"r")
        var = nffile.variables[varname]
        varvals = var.getValue()
        nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"a",mmap=False)
        var = nffile.variables[varname]
        var[:] = varvals
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"a")
        var = nffile.variables[varname]
        var.assignValue(varvals)
        nffile.close()
    ierr = 0
    return ierr

#======================================================================

UQdir = os.getcwd()

parm_names=[]
parm_indices=[]
parm_values=[]

myinput = open(options.parm_list, 'r')
lnum=0

username = getpass.getuser()
if (options.machine == 'cades' and options.runroot == ''):
    options.runroot = '/lustre/pfs1/cades-ccsi/scratch/'+username
elif (options.runroot == ''):
    options.runroot = '../../run'


#get parameter names and PFT information
casenames = []
if ('20TR' in options.casename or '1850' in options.casename):
    casenames.append(options.casename)
    isfullrun = False
else:
    casenames.append(options.casename+'_I1850CLM45CBCN_ad_spinup')
    casenames.append(options.casename+'_I1850CLM45CBCN')
    casenames.append(options.casename+'_I20TRCLM45CBCN')
    #for now, hard-code the number of years for ad_spinup and final spinup
    nyears_ad_spinup  = 250
    nyears_final_spinup = 250
    isfullrun = True

for s in myinput:
    pdata = s.split()
    #print pdata
    parm_names.append(pdata[0])
    parm_indices.append(int(pdata[1]))
myinput.close()

#get parameter values
myinput = open(options.parm_data, 'r')
for s in myinput:    
    parm_values.append(float(s))
myinput.close()

n_parameters = len(parm_names)
gst=str(100000+int(options.ensnum))

#create ensemble directories from original case(s)
isfirstcase = True
workdir = os.path.abspath('.')
for casename in casenames: 
    orig_dir = str(os.path.abspath(options.runroot)+'/'+casename+'/run')
    ens_dir  = os.path.abspath(options.runroot)+'/UQ/'+casename+'/g'+gst[1:]
    os.system('mkdir -p '+options.runroot+'/UQ/'+casename+'/g'+gst[1:]+'/timing/checkpoints')
    os.system('rm '+ens_dir+'/*.nc')
    os.system('cp '+orig_dir+'/*_in* '+ens_dir)
    os.system('cp '+orig_dir+'/*nml '+ens_dir)
    os.system('cp '+orig_dir+'/*stream* '+ens_dir)
    os.system('cp '+orig_dir+'/domain*.nc '+ens_dir)
    os.system('cp '+orig_dir+'/surf*.nc '+ens_dir)
    #os.system('cp '+orig_dir+'/*.r.*.nc '+ens_dir)
    os.system('cp '+orig_dir+'/*.rc '+ens_dir)
    os.system('cp '+orig_dir+'/*para*.nc '+ens_dir)
    os.system('cp '+orig_dir+'/*initial* '+ens_dir)
    os.system('cp '+orig_dir+'/*pftdyn* '+ens_dir)
    username = getpass.getuser()
    if (isfullrun):     #Full spinup simulation
        inifile=''
        if ('1850' in casename and 'ad_spinup' not in casename):
            yst = str(10000+nyears_ad_spinup+1)
            inifile = ens_dir_last+'/'+casename_last+'.clm2.r.'+yst[1:]+'-01-01-00000.nc'
        if ('20TR' in casename):
            yst = str(10000+nyears_final_spinup+1)
            inifile = ens_dir_last+'/'+casename_last+'.clm2.r.'+yst[1:]+'-01-01-00000.nc'
    else:                #Trasient case only
        inifile = ens_dir+'/'+casename+'.clm2.r.1974-01-01-00000.nc'
    casename_last = casename
    ens_dir_last = ens_dir


    #loop through all filenames, change directories in namelists, change parameter values
    for f in os.listdir(ens_dir):
        if (os.path.isfile(ens_dir+'/'+f) and (f[-2:] == 'in' or f[-3:] == 'nml' or 'streams' in f)):
            myinput=open(ens_dir+'/'+f)
            myoutput=open(ens_dir+'/'+f+'.tmp','w')
            for s in myinput:
                if ('paramfile' in s):
                    est = str(100000+int(options.ensnum))
                    os.system('cp '+ens_dir+'/clm_param*  '+ens_dir+'/clm_params_'+est[1:]+'.nc')
                    myoutput.write(" paramfile = './clm_params_"+est[1:]+".nc'\n")
                    #Hard-coded parameter file
                    pftfile = ens_dir+'/clm_params_'+est[1:]+'.nc'
                    pnum = 0
                    for p in parm_names:
                        #if (pnum == 0):
                        #    stem_leaf = getvar(pftfile, 'stem_leaf')
                        #    stem_leaf[2:5]=-1
                        #    ierr = putvar(pftfile, 'stem_leaf', stem_leaf)
                        param = getvar(pftfile, p)
                        if (parm_indices[pnum] > 0):
                            param[parm_indices[pnum]-1] = parm_values[pnum]
                        elif (parm_indices[pnum] == 0):
                            param = parm_values[pnum]
                        else:
                            param[:] = parm_values[pnum]
                        ierr = putvar(pftfile, p, param)
                        pnum = pnum+1
                #elif ('logfile =' in s):
                #    myoutput.write(s.replace('`date +%y%m%d-%H%M%S`',timestr))
                else:
                    myoutput.write(s.replace(orig_dir,ens_dir))
            myoutput.close()
            myinput.close()
            os.system(' mv '+ens_dir+'/'+f+'.tmp '+ens_dir+'/'+f)

    os.chdir(ens_dir)
    if (isfirstcase):
        exedir = os.path.abspath(orig_dir+'/../bld/')
    if (options.norun == False):
        if os.path.isfile(exedir+'/acme.exe'):
           os.system(exedir+'/acme.exe > acme_log.txt')
        elif os.path.isfile(exedir+'/e3sm.exe'):
           os.system(exedir+'/e3sm.exe > e3sm_log.txt')
        elif os.path.isfile(exedir+'/cesm.exe'):
           os.system(exedir+'/cesm.exe > cesm_log.txt')

    isfirstcase=False

#---------  code to post-process ensebmle member and cacluate total normalized SSE ----------
sse=0
myoutput = open('myoutput_sse.txt','w')
myind = 0
for p in parm_names:
        myoutput.write(str(parm_names[myind])+' '+str(parm_indices[myind])+' '+str(parm_values[myind])+'\n')
	myind = myind+1

for filename in os.listdir(UQdir+'/'+options.constraints):
  if (not os.path.isdir(filename)):
    myinput = open(UQdir+'/'+options.constraints+'/'+filename,'r')
    myvarname = filename.split('.')[0]  #model variable is filename
    #code to deal with special variables and/or aggregation
    #-------------
    lnum = 0
    year = 0
    for s in myinput:
        if (lnum == 0):
            header = s.split()
        else:
            hnum = 0
            PFT=-1      #default:  don't use PFT-specific info 
                        #  if specified, use h1 file (PFT-specific)
            doy=-1      #default:  annual average
            month=-1    #default:  don't use monthly data
            depth=-1
            unc = -999
            for h in header:
                if (h.lower() == 'year'):
                    year_info = s.split()[hnum]
                    if ('-' in year_info):
                        year_first = int(year_info.split('-')[0])
                        year_last  = int(year_info.split('-')[1])
                    else:
                        year_first = int(year_info)
                        year_last  = year_first
                if (h.lower() == 'doy'):
                    doy = int(s.split()[hnum])
                if (h.lower() == 'month'):
                    month = int(s.split()[hnum])
                if (h.lower() == 'pft'):
                    PFT = int(s.split()[hnum])
                if (h.lower() == 'value'):
                    value = float(s.split()[hnum])
                if (h.lower() == 'depth'):
                    depth = float(s.split()[hnum])
                if ('unc' in h.lower()):
                    unc   = float(s.split()[hnum])
                hnum = hnum+1
            #get the relevant variable/dataset
            #Assumes annual file(s) with daily output
            if (year <= 2016):
                if (PFT == -1):
                    filetype = 'h0'
                else:
                    filetype = 'h1'
                file_list = []
                for y in range(year_first,year_last+1):
                    yst_temp = str(10000+y)[1:]
                    file_list.append(casename+'.clm2.'+filetype+'.'+ \
                                    yst_temp+'-01-01-00000.nc')
                    
                #post processing of model output with nco to match constraining variables 
                for f in file_list:
                    if (myvarname == 'STEMC'):
                        os.system('ncap -s "STEMC=DEADSTEMC+LIVESTEMC" '+f+' '+f+'.tmp')
                        os.system('mv '+f+'.tmp '+f)
                        if (myvarname == 'AGBIOMASS'):
                            os.system('ncap -s "AGBIOMASS=DEADSTEMC+LIVESTEMC+LEAFC" '+f+' '+f+'.tmp')
                            os.system('mv '+f+'.tmp '+f)
                isfirstfile = True
                for f in file_list:
                    if (isfirstfile):
                        myvals = getvar(myfile, myvarname)/(year_last-year_first+1)
                        isfirstfile = False
                    else:
                        myvals = myvals + getvar(myfile, myvarname)/(year_last-year_first+1)

            if (doy > 0 and value > -900):
                if (PFT > 0):
                    #PFT-specific constraints 
                    model_val = myvals[doy,PFT-1]
                    if (unc < 0):
    	              unc = value*0.25 #default uncertainty set to 25%
                    sse = sse + ((model_val-value) /unc)**2
                    myoutput.write(str(myvarname)+' '+yst+' '+str(doy)+' '+str(PFT)+' '+ \
                                       str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
                elif (depth > 0):
                    #depth-specific constraint (column level only)
                    layers = [0,1.8,4.5,9.1,16.6,28.9,49.3,82.9,138.3,229.6,343.3]
                    for l in range(0,10):
                        if (depth >= layers[l] and depth < layers[l+1]):
                            thislayer = l
                            model_val = myvals[doy,thislayer,0]   
                            sse = sse + ((model_val-value) / unc )**2        
                            myoutput.write(str(myvarname)+' '+yst+' '+str(doy)+' '+str(depth)+' '+ \
                                               str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
                else:
                    #Daily column-level constraint, no depth/pft information (daily)
                    model_val = myvals[doy,0]
                    sse = sse + ((model_val-value) / unc )**2        
                    myoutput.write(str(myvarname)+' '+yst+' '+str(doy)+' '+str(PFT)+' '+ \
                                       str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
            elif (value > -900):
                if (PFT > 0):
                    #Annual constraints by PFT.  Assumes unit conversion from gC/m2/s to gC/m2/yr 
                    model_val = sum(myvals[0:,PFT-1])
                    if (myvarname == 'NPP' or myvarname == 'GPP' or myvarname == 'NEP' or myvarname == 'NEE'):
                        model_val = model_val * 24.0 * 3600  #convert to gC/m2/year
                    else:
                        model_val = model_val / 365.0        #mean value

                    sse = sse + ((model_val-value) / unc )**2      
                    myoutput.write(myvarname+' '+yst+' '+str(doy)+' '+str(PFT)+' '+ \
                                       str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
                else:
                    #Annual constraints (column level).  Assumes unit conversion from gC/m2/s to gC/m2/yr
                    model_val = sum(myvals[0:,0])
                    if (myvarname == 'NPP' or myvarname == 'GPP' or myvarname == 'NEP' or myvarname == 'NEE'):
                        model_val = model_val * 24.0 * 3600  #convert to gC/m2/year
                    else:
                        model_val = model_val / 365.0        #mean value
                    sse = sse + ((model_val-value) / unc )**2
                    myoutput.write(myvarname+' '+year_info+' '+str(doy)+' '+str(PFT)+' '+ \
                                       str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
        lnum = lnum+1
myoutput.close()

myoutput = open(workdir+'/qpso_ssedata/mysse_'+gst[1:]+'.txt','w')
myoutput.write(str(sse))
myoutput.close()



