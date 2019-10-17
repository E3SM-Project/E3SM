#!/usr/bin/env python

import netcdf_functions as nffun
import os, sys, csv, time, math, numpy
import re, subprocess
from optparse import OptionParser
#from Numeric import *


#pointCLM.py does the following:
#
#  1. Call routines to create surface and domain data (makepointdata.py)
#  2. Use create_newcase to build the new case with specified options
#  3. Set point and case-specific namelist options
#  4. configure case
#  5. build (compile) ACME with clean_build first if requested
#  6. apply user-specified PBS and submit information
#  7. submit single run or parameter ensemble job to PBS queue.
#
#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--diags", dest="diags", default=False, \
                 action="store_true", help="Write special outputs for diagnostics")
parser.add_option("--debug", dest="debug", default=False, \
                 action="store_true", help='Use debug queue and options')
parser.add_option("--runroot", dest="runroot", default="", \
                  help="Directory where the run would be created")
parser.add_option("--exeroot", dest="exeroot", default="", \
	         help="Location of executable")
parser.add_option("--lat_bounds", dest="lat_bounds", default='-999,-999', \
                  help = 'latitude range for regional run')
parser.add_option("--lon_bounds", dest="lon_bounds", default='-999,-999', \
                  help = 'longitude range for regional run')
parser.add_option("--humhol", dest="humhol", default=False, \
                  help = 'Use hummock/hollow microtopography', action="store_true")
parser.add_option("--marsh", dest="marsh", default=False, \
                  help = 'Use marsh hydrology/elevation', action="store_true")
parser.add_option("--mask", dest="mymask", default='', \
                  help = 'Mask file to use (regional only)')
parser.add_option("--ilambvars", dest="ilambvars", default=False, \
                 action="store_true", help="Write special outputs for diagnostics")
parser.add_option("--res", dest="res", default="CLM_USRDAT", \
                      help='Resoultion for global simulation')
parser.add_option("--point_list", dest="point_list", default='', \
                  help = 'File containing list of points to run')
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--coldstart", dest="coldstart", default=False, \
                  help = "set cold start (mutually exclusive w/finidat)", \
                  action="store_true")
parser.add_option("--compset", dest="compset", default='I1850CNPRDCTCBC', \
                  help = "component set to use (required)\n"
                         "Currently supports ONLY *CLM45(CN) compsets")
parser.add_option("--cruncep", dest="cruncep", default=False, \
                  help = "use cru-ncep data", action="store_true")
parser.add_option("--cplhist", dest="cplhist", default=False, \
                  help= "use CPLHIST forcing", action="store_true")
parser.add_option("--gswp3", dest="gswp3", default=False, \
                  help= "use GSWP3 forcing", action="store_true")
parser.add_option("--livneh", dest="livneh", default=False, \
                  action="store_true", help = "Livneh correction to CRU precip (CONUS only)")
parser.add_option("--machine", dest="machine", default = 'oic2', \
                  help = "machine to use (default = oic2)\n")
parser.add_option("--compiler", dest="compiler", default='gnu', \
	          help = "compiler to use (pgi, gnu)")
parser.add_option("--mpilib", dest="mpilib", default="mpi-serial", \
                      help = "mpi library (openmpi*, mpich, ibm, mpi-serial)")
parser.add_option("--ad_spinup", action="store_true", \
                  dest="ad_spinup", default=False, \
                  help = 'Run accelerated decomposition spinup')
parser.add_option("--exit_spinup", action="store_true", \
                  dest="exit_spinup", default=False, \
                  help = 'Run exit spinup (CLM 4.0 only)')
parser.add_option("--csmdir", dest="csmdir", default='/../../../../..', \
                  help = "base CESM directory (default = ../../..)")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../../../../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--finidat_case", dest="finidat_case", default='', \
                  help = "case containing initial data file to use" \
                  +" (should be in your run directory)")
parser.add_option("--finidat", dest="finidat", default='', \
                  help = "initial data file to use" \
                  +" (should be in your run directory)")
parser.add_option("--finidat_year", dest="finidat_year", default=-1, \
                  help = "model year of initial data file (default is" \
                  +" last available)")
parser.add_option("--run_units", dest="run_units", default='nyears', \
                  help = "run length units (ndays, nyears)")
parser.add_option("--run_n", dest="run_n", default=50, \
                  help = "run length (in run units)")
parser.add_option("--run_startyear", dest="run_startyear",default=-1, \
                      help='Starting year for model output')
parser.add_option("--rmold", dest="rmold", default=False, action="store_true", \
                  help = 'Remove old case directory with same name' \
                  +" before proceeding")
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--parm_file", dest="parm_file", default='',
                  help = 'file for parameter modifications')
parser.add_option("--parm_vals", dest="parm_vals", default="", \
                  help = 'User specified parameter values')
parser.add_option("--parm_file_P", dest="parm_file_P", default='',
                  help = 'file for P parameter modifications')
parser.add_option("--hist_mfilt", dest="hist_mfilt", default=-1, \
                  help = 'number of output timesteps per file')
parser.add_option("--hist_nhtfrq", dest="hist_nhtfrq", default=-999, \
                  help = 'output file timestep')
parser.add_option("--hist_vars", dest="hist_vars", default='', \
                  help = 'use hist_vars file')
#parser.add_option("--queue", dest="queue", default='essg08q', \
#                  help = 'PBS submission queue')
parser.add_option("--clean_config", dest="clean_config", default=False, \
                  help = 'Run cesm_setup -clean script')
parser.add_option("--clean_build", dest="clean_build", default=False, \
                  help = 'Perform clean build before building', \
                  action="store_true")
parser.add_option("--no_config", dest="no_config", default=False, \
                  help = 'do NOT configure case', action="store_true")
parser.add_option("--no_build", dest="no_build", default=False, \
                  help = 'do NOT build CESM', action="store_true")
parser.add_option("--no_submit", dest="no_submit", default=False, \
                  help = 'do NOT submit CESM to queue', action="store_true")
parser.add_option("--align_year", dest="align_year", default=-999, \
                  help = 'Alignment year (transient run only)')
parser.add_option("--np", dest="np", default=1, \
                  help = 'number of processors')
parser.add_option("--ninst", dest="ninst", default=1, \
                  help = 'number of land model instances')
parser.add_option("--ng", dest="ng", default=64, \
                  help = 'number of groups to run in ensmble mode')
parser.add_option("--tstep", dest="tstep", default=0.5, \
                  help = 'CLM timestep (hours)')
parser.add_option("--co2_file", dest="co2_file", default="fco2_datm_1765-2007_c100614.nc", \
                  help = 'CLM timestep (hours)')
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=250, \
                  help = 'number of years to run ad_spinup')
parser.add_option("--metdir", dest="metdir", default="none", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--nopointdata", action="store_true", \
                  dest="nopointdata", help="Do NOT make point data (use data already created)", \
                  default=False)
#parser.add_option("--cleanlogs",dest="cleanlogs", help=\
#                   "Removes temporary and log files that are created",\
#                   default=False,action="store_true")
parser.add_option("--nofire", action="store_true", dest="nofire", default=False, \
                    help="To turn off wildfires")
parser.add_option("--nopftdyn", action="store_true", dest="nopftdyn", \
                      default = False, help='Do not use dynamic PFT file')
parser.add_option("--harvmod", action="store_true", dest="harvmod", \
                      default=False, help = "Turn on harvest modificaton" \
                      "All harvest is performed in first timestep")
parser.add_option("--no_dynroot", dest="no_dynroot", default=False, \
                  help = 'Turn off dynamic root distribution', action="store_true")
parser.add_option("--bulk_denitrif", dest="bulk_denitrif", default=False, \
                  help = 'To turn off BGC nitrification-denitrification', action="store_true")
parser.add_option("--vertsoilc", dest="vsoilc", default=False, \
                  help = 'To turn on CN with multiple soil layers, excluding CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me', action="store_true")
parser.add_option("--C13", dest="C13", default=False, \
                  help = 'Switch to turn on C13', action="store_true")
parser.add_option("--C14", dest="C14", default=False, \
                  help = 'Use C14 as C13 (no decay)', action="store_true")
parser.add_option("--branch", dest="branch", default=False, \
		  help = 'Switch for branch run', action="store_true")
parser.add_option("--makemetdata", dest="makemet", default=False, \
		  help = 'Generate meteorology', action="store_true")
parser.add_option("--surfdata_grid", dest="surfdata_grid", default=False, \
                  help = 'Use gridded surface data instead of site data', action="store_true")
parser.add_option("--include_nonveg", dest="include_nonveg", default=False, \
                  help = 'Include non-vegetated columns/Landunits in surface data')
parser.add_option("--trans2", dest="trans2", default=False, action="store_true", \
                  help = 'Tranisnent phase 2 (1901-2010) - CRUNCEP only')
parser.add_option("--spinup_vars", dest="spinup_vars", default=False, \
                  help = 'Limit output vars in spinup runs', action="store_true")
parser.add_option("--trans_varlist", dest = "trans_varlist", default='', help = "Transient outputs")
parser.add_option("--c_only", dest="c_only", default=False, \
                 help="Carbon only (supplemental P and N)", action="store_true")
parser.add_option("--cn_only", dest="cn_only", default=False, \
                  help = 'Carbon/Nitrogen only (supplemental P)', action="store_true")
parser.add_option("--cp_only", dest="cp_only", default=False, \
                  help = 'Carbon/Phosphorus only (supplemental N)', action = "store_true")
parser.add_option("--ensemble_file", dest="ensemble_file", default='', \
                  help = 'Parameter sample file to generate ensemble')
parser.add_option("--mc_ensemble", dest="mc_ensemble", default=-1, \
                  help = 'Monte Carlo ensemble (argument is # of simulations)')
parser.add_option("--ensemble_nocopy", dest="ensemble_nocopy", default=False, \
                  help = 'Do not copy files to ensemble directories', action="store_true")
#Changed by Ming for mesabi
parser.add_option("--archiveroot", dest="archiveroot", default='', \
                  help = "archive root directory only for mesabi")
#Added by Kirk to include the modified parameter file
parser.add_option("--mod_parm_file", dest="mod_parm_file", default='', \
                  help = "adding the path to the modified parameter file")
parser.add_option("--mod_parm_file_P", dest="mod_parm_file_P", default='', \
                  help = "adding the path to the modified parameter file")
parser.add_option("--parm_list", dest="parm_list", default='parm_list', \
                  help = 'File containing list of parameters to vary')
parser.add_option("--postproc_file", dest="postproc_file", default="", \
                  help = 'File for ensemble post processing')
parser.add_option("--walltime", dest="walltime", default=6, \
                  help = "desired walltime for each job (hours)")
(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------

#machine info:  cores per node
ppn=1
if ('titan' in options.machine):
    ppn=16
    if (int(options.walltime) > 2 and int(options.ng) < 2048):
        print 'Requested walltime too long'
        print 'Setting to 2 hours.'
        options.walltime=2
elif ('metis' in options.machine):
    ppn=16
elif ('oic2' in options.machine):
    ppn=8
elif ('oic5' in options.machine or 'cori-haswell' in options.machine or 'eos' in options.machine \
      or 'cades' in options.machine):
    ppn=32
elif ('cori-knl' in options.machine):
    ppn=64
elif ('edison' in options.machine):
    ppn=24

PTCLMdir = os.getcwd()

if (options.hist_vars != ''):
    hist_vars = os.path.abspath(options.hist_vars)

#check for valid csm directory
if (os.path.exists(options.csmdir) == False):
    print('Error:  invalid CESM root directory')
    sys.exit()
else:
    csmdir=os.path.abspath(PTCLMdir+'/'+options.csmdir)
    scriptsdir = csmdir+'/cime/scripts'
#case directory
if (options.caseroot == '' or (os.path.exists(options.caseroot) == False)):
    caseroot = csmdir+'/cime/scripts'
else:
    caseroot = os.path.abspath(options.caseroot)

#case run root directory
if (options.runroot == '' or (os.path.exists(options.runroot) == False)):
    runroot = csmdir+'/run'
else:
    runroot = os.path.abspath(options.runroot)

#check for valid input data directory
print(options.ccsm_input)
if (options.ccsm_input == '' or (os.path.exists(options.ccsm_input) \
                                 == False)):
    print('Error:  invalid input data directory')
    sys.exit()
else:
    options.ccsm_input = os.path.abspath(options.ccsm_input)

compset = options.compset
isglobal = False
if (options.site == ''):
    isglobal = True
    options.site=options.res

if ('CBCN' in compset or 'ICB' in compset or 'CLM45CB' in compset):
    cpl_bypass = True
else:
    cpl_bypass = False

surfdir = 'surfdata_map'
#pftphys_stamp = '_c160711' #'c160711_root'
#pftphys_stamp = '_c160711_test170303'
pftphys_stamp = '_c170808'
if ('ECA' in compset):
    pftphys_stamp = '.c160709'
#CNPstamp = 'c131108'
CNPstamp = 'c170306'

#check consistency of options
if ('20TR' in compset):
    #ignore spinup option if transient compset
    if (options.ad_spinup or options.exit_spinup):
      print('Spinup options not available for transient compset.')
      sys.exit()
    #finidat is required for transient compset
    if (options.finidat_case == '' and options.finidat == ''):
        print('Error:  must provide initial data file for I20TR compsets')
        sys.exit()

#get full path of finidat file
finidat=''
finidat_year=int(options.finidat_year)

if ('CN' in compset or 'ECA' in compset):
  mybgc = 'CN'
elif ('ED' in compset):
  mybgc = 'ED'
else:
  mybgc = 'none'

if (options.exit_spinup):
    if (options.mycaseid != ''):
        finidat = options.mycaseid+'_'+options.site+'_I1850'+mybgc+'_ad_spinup'
    else:
        finidat = options.site+'_I1850'+mybgc+'_ad_spinup'
    finidat_year = int(options.ny_ad)+1

if (options.finidat == ''  and options.finidat_case == ''):  #not user-defined
    if (options.coldstart==False and compset == "I1850CLM45"+mybgc and options.ad_spinup == False):
        if (options.mycaseid != ''):
            options.finidat_case = options.mycaseid+'_'+options.site+ \
                '_I1850CLM45'+mybgc+'_ad_spinup'
        else:
            options.finidat_case = options.site+'_I1850CLM45'+mybgc+'_ad_spinup'

        if (options.finidat_year == -1):
            finidat_year = int(options.ny_ad)+1
    
    if (compset == "I20TRCLM45"+mybgc or compset == "I20TRCRUCLM45"+mybgc):
        if (options.mycaseid != ''):
            options.finidat_case = options.mycaseid+'_'+options.site+ \
                '_I1850CLM45'+mybgc
        else:
            options.finidat_case = options.site + '_I1850CLM45'+mybgc
            
        if (options.finidat_year == -1):
            finidat_year=1850
    if (compset == "I20TR"+mybgc):
        if (options.mycaseid != ''):
            options.finidat_case = options.mycaseid+'_'+options.site+ \
                '_I1850'+mybgc
        else:
            options.finidat_case = options.site + '_I1850'+mybgc
        if (options.finidat_year == -1):
            finidat_year = 1850

        #finidat is required for transient compset
            if (os.path.exists(runroot+'/'+options.finidat_case) == False):
                print('Error:  must provide initial data file for I20TRCLM45CN/BGC compset OR '+ \
                          runroot+'/'+options.finidat_case+' existed as refcase')
                sys.exit


if (options.finidat_case != ''):
    finidat_yst = str(10000+finidat_year)
    finidat = runroot+'/'+options.finidat_case+'/run/'+ \
              options.finidat_case+'.clm2.r.'+finidat_yst[1:]+ \
              '-01-01-00000.nc'

if (options.finidat != ''):
    finidat = options.finidat
    finidat_year = int(finidat[-19:-15])
    finidat_yst = str(10000+finidat_year)    

#construct default casename
casename    = options.site+"_"+compset
if (options.mycaseid != ""):
    casename = options.mycaseid+'_'+casename

#CRU-NCEP 2 transient phases
if ('CRU' in compset or options.cruncep or options.gswp3 or options.cplhist):
    use_reanalysis = True
else:
    use_reanalysis = False
if ('20TR' in compset and use_reanalysis and (not cpl_bypass)):
    if options.trans2:
        casename = casename+'_phase2'
    else:
        casename = casename+'_phase1'
if (options.ad_spinup):
    casename = casename+'_ad_spinup'
if (options.exit_spinup):
    casename = casename+'_exit_spinup'

PTCLMfiledir = options.ccsm_input+'/lnd/clm2/PTCLM'

if (caseroot != "./"):
    casedir=caseroot+"/"+casename
else:
    casedir=casename

#Check for existing case directory
if (os.path.exists(casedir)):
    
    print('Warning:  Case directory exists')
    if (options.rmold):
        print('--rmold specified.  Removing old case ')
        os.system('rm -rf '+casedir)
    else:
        var = raw_input('proceed (p), remove old (r), or exit (x)? ')
        if var[0] == 'r':
            os.system('rm -rf '+casedir)
        if var[0] == 'x':
            sys.exit()    
print("CASE directory is: "+casedir+"\n")

#Construct case build and run directory
if (options.exeroot == '' or (os.path.exists(options.exeroot) == False)):
    exeroot = runroot+'/'+casename+'/bld'
    #if ('titan' in options.machine or 'eos' in options.machine):
    #    exeroot = os.path.abspath(os.environ['HOME']+ \
   # 	    '/acme_scratch/pointclm/'+casename+'/bld')
else:
    exeroot=options.exeroot
print("CASE exeroot is: "+exeroot+"\n")
rundir=runroot+'/'+casename+'/run'
print("CASE rundir is: "+rundir+"\n")
if (options.rmold):
    if (options.no_build == False):
        print('Removing build directory: '+exeroot)
        os.system('rm -rf '+exeroot)
    print('Removing run directory: '+rundir)
    os.system('rm -rf '+rundir)

#------Make domain, surface data and pftdyn files ------------------
mysimyr=1850
if ('1850' not in compset and '20TR' not in compset):
    mysimyr=2000
print compset, mysimyr

if (options.nopointdata == False):
    ptcmd = 'python makepointdata.py --caseroot '+caseroot+ \
        ' --csmdir '+csmdir+' --ccsm_input '+options.ccsm_input+ \
        ' --lat_bounds '+options.lat_bounds+' --lon_bounds '+ \
        options.lon_bounds+' --mysimyr '+str(mysimyr)
    if (options.metdir != 'none'):
        ptcmd = ptcmd + ' --metdir '+options.metdir
    if (options.makemet):
        ptcmd = ptcmd + ' --makemetdata'
    if (options.surfdata_grid):
        ptcmd = ptcmd + ' --surfdata_grid'
    if (options.include_nonveg):
        ptcmd = ptcmd + ' --include_nonveg'
    if (options.nopftdyn):
        ptcmd = ptcmd + ' --nopftdyn'
    if (options.mymask != ''):
        ptcmd = ptcmd + ' --mask '+options.mymask
    if (isglobal):
        ptcmd = ptcmd + ' --res '+options.res
        if (options.point_list != ''):
            ptcmd = ptcmd+' --point_list '+options.point_list

    else:
        ptcmd = ptcmd + ' --site '+options.site+' --sitegroup '+options.sitegroup       

    print('Creating surface and domain data')
    if (options.machine == 'eos' or options.machine == 'titan'):
        os.system('rm temp/*.nc')
        print('Note:  NCO operations are slow on eos and titan.')
        print('Submitting PBS script to make surface and domain data on rhea')
        pbs_rhea=open('makepointdata_rhea.pbs','w')
        pbs_rhea.write('#PBS -l walltime=00:30:00\n')
        pbs_rhea.write('#PBS -l nodes=1\n')
        pbs_rhea.write('#PBS -A cli112\n')
        pbs_rhea.write('#PBS -q rhea\n')
        pbs_rhea.write('#PBS -l gres=atlas1%atlas2\n\n')
        pbs_rhea.write(' module load nco\n')
        pbs_rhea.write(' cd '+csmdir+'/components/clm/tools/clm4_5/pointclm/\n')
        pbs_rhea.write(' rm temp/*.nc\n')
        pbs_rhea.write(' module unload PE-intel\n')
        pbs_rhea.write(' module load PE-gnu\n')
        pbs_rhea.write(' module load python\n')
        pbs_rhea.write(' module load python_numpy\n')
        pbs_rhea.write(' module load python_scipy\n')
        pbs_rhea.write(ptcmd+'\n')
        pbs_rhea.close()
        os.system('qsub makepointdata_rhea.pbs')
        n_nc_files = 3
        if (options.nopftdyn):
            n_nc_files = 2
        n=0
        while (n < n_nc_files):
            #Wait until files have been generated on rhea to proceed
            list_dir = os.listdir('./temp')
            n=0
            for file in list_dir:
                if file.endswith('.nc'):
                    n=n+1
            os.system('sleep 10')
        #Clean up
        os.system('rm makepointdata_rhea*') 
    else:
        os.system(ptcmd)


#get site year information
sitedatadir = os.path.abspath(PTCLMfiledir)
os.chdir(sitedatadir)
if (isglobal == False):
    AFdatareader = csv.reader(open(options.sitegroup+'_sitedata.txt',"rb"))
    for row in AFdatareader:
        if row[0] == options.site:
            if (use_reanalysis):
                if ('CN' in compset or 'BGC' in compset):
                    if (options.trans2):
                        startyear = 1921
                        endyear   = int(row[7])
                    else:
                        startyear = 1901
                        endyear   = 1920
                else:
                    startyear = int(row[6]) #1901
                    endyear   = int(row[7])
            else:
                startyear=int(row[6])
                endyear=int(row[7])
            alignyear = int(row[8])
            if (options.diags):
                timezone = int(row[9])
            if ('US-SPR' in options.site):
              numxpts=2
            else:
              numxpts=1
            numypts=1
else:
    if (use_reanalysis):
        startyear=1901
        if ('20TR' in compset):
            endyear = 2010
            if (options.trans2):
                startyear = 1921
        else:
            endyear   = 1920
    else:    #Global default to Qian
        startyear=1948
        if ('20TR' in compset):
            endyear = 2004
            if (options.trans2):
                startyear = 1973
        else:
            endyear   = 1972

if (isglobal == False):
    ptstr = str(numxpts)+'x'+str(numypts)+'pt'
os.chdir(csmdir+'/cime/scripts')

#parameter (pft-phys) modifications if desired
tmpdir = csmdir+'/components/clm/tools/clm4_5/pointclm/temp'
os.system('mkdir -p '+tmpdir)
if (options.mod_parm_file != ''):
    os.system('nccopy -3 '+options.mod_parm_file+' '+tmpdir+'/clm_params.nc')
else:
    os.system('nccopy -3 '+options.ccsm_input+'/lnd/clm2/paramdata/clm_params'+pftphys_stamp+'.nc ' \
              +tmpdir+'/clm_params.nc')
os.system('chmod u+w ' +tmpdir+'/clm_params.nc')
if (options.parm_file != ''):
    pftfile = tmpdir+'/clm_params.nc'
    if ('/' not in options.parm_file):
       #assume in pointclm directory
       input  = open(PTCLMdir+'/'+options.parm_file)
    else:   #assume full path given
       input   = open(os.path.abspath(options.parm_file))
    for s in input:
        if s[0:1] != '#':
            values = s.split()
            thisvar = nffun.getvar(pftfile, values[0])
            if (len(values) == 2):
                if isinstance(thisvar,list):
                    thisvar[:] = float(values[1])
                else:
                    thisvar = float(values[1])  #scalar
            elif (len(values) == 3):
                if (float(values[1]) > 0):
                    thisvar[int(values[1])] = float(values[2])
                else:
                    if isinstance(thisvar,list):
                        thisvar[:] = float(values[2])
                    else:
                        thisvar = float(values[2])
            ierr = nffun.putvar(pftfile, values[0], thisvar)
    input.close()
if (options.parm_vals != ''):
    pftfile = tmpdir+'/clm_params.nc'
    parms = options.parm_vals.split('/')
    nparms = len(parms)
    for n in range(0,nparms):
	parm_data=parms[n].split(',')
        thisvar = nffun.getvar(pftfile, parm_data[0])
        if (len(parm_data) == 2):
	   thisvar[:] = float(parm_data[1])
        elif (len(parm_data) == 3): 
           if (float(parm_data[1]) >= 0):
               thisvar[int(parm_data[1])] = float(parm_data[2])
           else: 
               thisvar[:] = float(parm_data[2])
        ierr =  nffun.putvar(pftfile, parm_data[0], thisvar)
           
#parameter (soil order dependent) modifications if desired    ::X.YANG 
if (options.mod_parm_file_P != ''):
    os.system('cp '+options.mod_parm_file_P+' '+tmpdir+'/CNP_parameters.nc')
else:
    os.system('cp '+options.ccsm_input+'/lnd/clm2/paramdata/CNP_parameters_'+CNPstamp+'.nc ' \
              +tmpdir+'/CNP_parameters.nc')
os.system('chmod u+w ' +tmpdir+'/CNP_parameters.nc')


if (options.parm_file_P != ''):
    soilorderfile = tmpdir+'/CNP_parameters.nc'
    if ('/' not in options.parm_file_P):
        #assume in pointclm directory
        input  = open(PTCLMdir+'/'+options.parm_file_P)
    else:   #assume full path given
        input   = open(os.path.abspath(options.parm_file_P))
        input   = open(os.path.abspath(options.parm_file_P))
    for s in input:
        if s[0:1] != '#':
            values = s.split()
            thisvar = nffun.getvar(soilorderfile, values[0])
            if (len(values) == 2):
                thisvar[:] = float(values[1])
            elif (len(values) == 3):
                if (float(values[1]) >= 0):
                    thisvar[int(values[1])] = float(values[2])
                else:
                    thisvar[:] = float(values[2])
            ierr = nffun.putvar(soilorderfile, values[0], thisvar)
    input.close()

#set number of run years for ad, exit spinup cases
if (options.ny_ad != options.run_n and options.ad_spinup):
    options.run_n = options.ny_ad
elif (options.exit_spinup):
    options.run_n = 1

#create new case
print ('./create_newcase --case '+casename+' --mach '+options.machine+' --compset '+ \
           options.compset+' --res '+options.res+' --compiler '+options.compiler+' --mpilib '+ \
           options.mpilib+' --walltime '+str(options.walltime)+':00:00')
os.system('./create_newcase --case '+casename+' --mach '+options.machine+' --compset '+ \
	   options.compset+' --res '+options.res+' --compiler '+options.compiler+' --mpilib '+ \
           options.mpilib+' > create_newcase.log'+' --walltime '+str(options.walltime)+ \
          ':00:00')
if (os.path.isdir(casename)):
    print(casename+' created.  See create_newcase.log for details')
    os.system('mv create_newcase.log '+casename)
else:
    print('failed to create case.  See create_newcase.log for details')

os.chdir(casedir)

os.system('./xmlchange SAVE_TIMING=FALSE')
#------------------ env_build.xml modifications ------------------------
os.system('./xmlchange -file env_build.xml -id EXEROOT -val '+exeroot)

#turn off ROF module
os.system('./xmlchange -file env_build.xml -id RTM_MODE -val NULL')

#clm 4_5 cn config options
clmcn_opts = "-phys clm4_5"

os.system('./xmlchange -file env_build.xml -id CLM_CONFIG_OPTS -val "'+ \
              clmcn_opts+'"')
print("CLM module options: " + clmcn_opts+"\n")
if (options.machine == 'userdefined'):
    os.system('./xmlchange -file env_build.xml -id COMPILER -val "' + \
                  options.compiler+'"')
    os.system('./xmlchange -file env_build.xml -id OS -val "' + \
                  'linux"')
    os.system('./xmlchange -file env_build.xml -id EXEROOT -val "'+runroot+'/'+casename+'/bld"')

#-------------- env_run.xml modifications -------------------------
#if (options.runroot != ''):
os.system('./xmlchange -file env_run.xml -id RUNDIR -val '+rundir)
os.system('./xmlchange -file env_run.xml -id DOUT_S -val TRUE')
os.system('./xmlchange -file env_run.xml -id DOUT_S_ROOT -val ' \
            +runroot+'/archive/'+casename)
if (options.ccsm_input != ''):
    os.system('./xmlchange -file env_run.xml -id DIN_LOC_ROOT -val ' \
                  +options.ccsm_input)
    
#define mask and resoultion
if (isglobal == False):
    os.system('./xmlchange -file env_run.xml -id CLM_USRDAT_NAME ' \
                  +' -val '+str(numxpts)+'x'+str(numypts)+'pt_'+options.site)
if (options.ad_spinup):
    os.system('./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS ' \
                  +'-append -val "-bgc_spinup on"')
if (options.run_startyear > -1):
    os.system('./xmlchange -file env_run.xml -id RUN_STARTDATE ' \
                  +' -val '+str(options.run_startyear)+'-01-01')
os.system('./xmlchange -file env_run.xml -id ATM_DOMAIN_PATH ' \
              +' -val "\${RUNDIR}"')
os.system('./xmlchange -file env_run.xml -id LND_DOMAIN_PATH ' \
              +' -val "\${RUNDIR}"')
os.system('./xmlchange -file env_run.xml -id ATM_DOMAIN_FILE ' \
              +' -val domain.nc')
os.system('./xmlchange -file env_run.xml -id LND_DOMAIN_FILE ' \
              +' -val domain.nc')

#turn off archiving
os.system('./xmlchange -file env_run.xml -id DOUT_S ' \
                  +' -val "FALSE"') 
#datm options
if (not cpl_bypass):
    if (use_reanalysis):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'DATM_MODE -val CLMCRUNCEP') 
    else:
        if (isglobal == False):
            os.system('./xmlchange -file env_run.xml -id ' \
                          +'DATM_MODE -val CLM1PT') 
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_START -val '+str(startyear))
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_END -val '+str(endyear))
    if (options.align_year == -999):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'DATM_CLMNCEP_YR_ALIGN -val '+str(1))
    else:
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_ALIGN -val '+str(options.align_year))


#Set input data directory
os.system('./xmlchange -file env_run.xml -id ' \
              +'DIN_LOC_ROOT -val '+options.ccsm_input)

#Change simulation timestep
if (options.tstep != 0.5):
    os.system('./xmlchange -file env_run.xml -id ' \
                      +'ATM_NCPL -val '+str(int(24/float(options.tstep))))

#Branch run options
if (options.branch or options.exit_spinup):
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'RUN_TYPE -val branch')
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'RUN_REFDATE -val '+finidat_yst[1:]+'-01-01')
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'RUN_REFCASE -val '+options.finidat_case)
else:
    if (('CN' in compset or 'BGC' in compset) and options.ad_spinup \
            == False and options.coldstart==False):
        os.system('./xmlchange -file env_run.xml -id RUN_REFDATE -val ' \
                      +finidat_yst[1:]+'-01-01')

    #adds capability to run with transient CO2
if ('20TR' in compset):
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'CCSM_BGC -val CO2A')
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'CLM_CO2_TYPE -val diagnostic')
    if (options.run_startyear == -1):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_STARTDATE -val 1850-01-01')
    
#no PIO on oic
#if ('oic' in options.machine or 'eos' in options.machine or 'edison' in options.machine):
os.system('./xmlchange -file env_run.xml -id PIO_TYPENAME -val netcdf')

os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val '+str(options.np))
if (options.np == 1):
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val 1')
if (int(options.np) > 1):
    os.system('./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE -val '+str(ppn))
    os.system('./xmlchange -file env_mach_pes.xml -id MAX_MPITASKS_PER_NODE -val '+str(ppn))

if (int(options.ninst) > 1):
    os.system('./xmlchange -file env_mach_pes.xml -id ' \
                  +'NINST_LND -val '+options.ninst)
    os.system('./xmlchange -file env_mach_pes.xml -id ' \
                  +'NTASKS_LND -val '+options.ninst)

os.system('./xmlchange -file env_run.xml -id ' \
              +'STOP_OPTION -val '+options.run_units)
os.system('./xmlchange -file env_run.xml -id ' \
              +'STOP_N -val '+str(options.run_n))


#--------------------------CESM setup ----------------------------------------

if (options.clean_config):
    os.system('./case.setup -clean')
    os.system('rm -f Macro')
    os.system('rm -f user-nl-*')

# Add options for FFLAGS to Macros file here 

#clm namelist modifications
for i in range(1,int(options.ninst)+1):
    if (int(options.ninst) == 1):
        output = open("user_nl_clm",'w')
    else:
        if (i < 10):
            output = open("user_nl_clm_000"+str(i),'w')
        elif (i < 100):
            output = open("user_nl_clm_00"+str(i),'w')
        elif (i < 1000):
            output = open("user_nl_clm_0"+str(i),'w')
    output.write('&clm_inparm\n')

    #history file options
    #outputs for SPRUCE MiP and Jiafu's diagnostics code:
    var_list_hourly = ['GPP', 'NEE', 'NEP', 'NPP', 'LEAFC_ALLOC', 'AGNPP', 'MR', \
            'CPOOL_TO_DEADSTEMC', 'LIVECROOTC_XFER_TO_LIVECROOTC', 'DEADCROOTC_XFER_TO_DEADCROOTC', \
            'CPOOL_TO_LIVECROOTC', 'CPOOL_TO_DEADCROOTC', 'FROOTC_ALLOC', 'AR', 'LEAF_MR', 'CPOOL_LEAF_GR',
            'TRANSFER_LEAF_GR', 'CPOOL_LEAF_STORAGE_GR', 'LIVESTEM_MR', 'CPOOL_LIVESTEM_GR', \
            'TRANSFER_LIVESTEM_GR', 'CPOOL_LIVESTEM_STORAGE_GR', 'CPOOL_DEADSTEM_GR', \
            'TRANSFER_DEADSTEM_GR', 'CPOOL_DEADSTEM_STORAGE_GR', 'LIVECROOT_MR', 'CPOOL_LIVECROOT_GR', \
            'TRANSFER_LIVECROOT_GR', 'CPOOL_LIVECROOT_STORAGE_GR', 'CPOOL_DEADCROOT_GR', 'TRANSFER_DEADCROOT_GR', 'CPOOL_DEADCROOT_STORAGE_GR', \
            'FROOT_MR', 'CPOOL_FROOT_GR', 'TRANSFER_FROOT_GR', 'CPOOL_FROOT_STORAGE_GR', 'FSH', 'EFLX_LH_TOT', \
            'Rnet', 'FCTR', 'FGEV', 'FCEV', 'SOILLIQ', 'QOVER', 'QDRAI', 'TOTVEGC', 'LEAFC', 'LIVESTEMC', 'DEADSTEMC', \
            'FROOTC', 'LIVECROOTC', 'DEADCROOTC', 'TG', 'TV', 'TSA', 'TSOI', 'DEADSTEMC_STORAGE', \
            'LIVESTEMC_STORAGE', 'DEADCROOTC_STORAGE', 'LIVECROOTC_STORAGE', 'CPOOL_TO_DEADSTEMC_STORAGE', \
            'CPOOL_TO_LIVESTEMC_STORAGE', 'CPOOL_TO_DEADCROOTC_STORAGE', 'CPOOL_TO_LIVECROOTC_STORAGE', \
            'ER', 'HR', 'FROOTC_STORAGE', 'LEAFC_STORAGE', 'LEAFC_XFER', 'FROOTC_XFER', 'LIVESTEMC_XFER', \
            'DEADSTEMC_XFER', 'LIVECROOTC_XFER', 'DEADCROOTC_XFER', 'SR', 'HR_vr', 'FIRA', 
            'FSA', 'FSDS', 'FLDS', 'TBOT', 'RAIN', 'SNOW', 'WIND', 'PBOT', 'QBOT', 'QVEGT', 'QVEGE', 'QSOIL', \
            'QFLX_SUB_SNOW', 'QFLX_DEW_GRND', 'QH2OSFC', 'H2OSOI', 'CPOOL_TO_LIVESTEMC', 'TOTLITC', \
            'TOTSOMC', 'ZWT', 'SNOWDP', 'TLAI','RH2M']
    var_list_daily = ['TOTLITC', 'TOTSOMC', 'CWDC', 'LITR1C_vr', 'LITR2C_vr', 'LITR3C_vr', 'SOIL1C_vr', \
                      'SOIL2C_vr', 'SOIL3C_vr', 'H2OSFC', 'ZWT', 'SNOWDP', 'TLAI', 'CPOOL','NPOOL','PPOOL', \
                      'FPI','FPI_P','FPG','FPG_P','FPI_vr','FPI_P_vr']
    var_list_pft = ['GPP', 'NPP', 'LEAFC_ALLOC', 'AGNPP', 'CPOOL_TO_DEADSTEMC', \
                    'LIVECROOTC_XFER_TO_LIVECROOTC', 'DEADCROOTC_XFER_TO_DEADCROOTC', \
                    'CPOOL_TO_LIVECROOTC', 'CPOOL_TO_DEADCROOTC', 'FROOTC_ALLOC', 'AR', 'MR', \
                    'LEAF_MR', 'CPOOL_LEAF_GR', 'TRANSFER_LEAF_GR', 'CPOOL_LEAF_STORAGE_GR', \
                    'LIVESTEM_MR', 'CPOOL_LIVESTEM_GR', 'TRANSFER_LIVESTEM_GR', \
                    'CPOOL_LIVESTEM_STORAGE_GR', 'CPOOL_DEADSTEM_GR', 'TRANSFER_DEADSTEM_GR', \
                    'CPOOL_DEADSTEM_STORAGE_GR', 'LIVECROOT_MR', 'CPOOL_LIVECROOT_GR', \
                    'TRANSFER_LIVECROOT_GR', 'CPOOL_LIVECROOT_STORAGE_GR', 'CPOOL_DEADCROOT_GR', \
                    'TRANSFER_DEADCROOT_GR', 'CPOOL_DEADCROOT_STORAGE_GR', 'FROOT_MR', \
                    'CPOOL_FROOT_GR', 'TRANSFER_FROOT_GR', 'CPOOL_FROOT_STORAGE_GR', 'FCTR', 'FCEV', \
                    'TOTVEGC', 'LEAFC', 'LIVESTEMC', 'DEADSTEMC', 'FROOTC', 'LIVECROOTC', \
                    'DEADCROOTC', 'DEADSTEMC_STORAGE', 'LIVESTEMC_STORAGE', 'DEADCROOTC_STORAGE', \
                    'LIVECROOTC_STORAGE', 'CPOOL_TO_DEADSTEMC_STORAGE', 'CPOOL_TO_LIVESTEMC_STORAGE', \
                    'CPOOL_TO_DEADCROOTC_STORAGE', 'CPOOL_TO_LIVECROOTC_STORAGE', \
                    'FROOTC_STORAGE', 'LEAFC_STORAGE', 'LEAFC_XFER', 'FROOTC_XFER', 'LIVESTEMC_XFER', \
                    'DEADSTEMC_XFER', 'LIVECROOTC_XFER', 'DEADCROOTC_XFER', 'TLAI', 'CPOOL_TO_LIVESTEMC']
    var_list_spinup = ['PPOOL', 'EFLX_LH_TOT', 'RETRANSN', 'PCO2', 'PBOT', 'NDEP_TO_SMINN', 'OCDEP', \
                       'BCDEP', 'COL_FIRE_CLOSS', 'HDM', 'LNFM', 'NEE', 'GPP', 'FPSN', 'AR', 'HR', \
                       'MR', 'GR', 'ER', 'NPP', 'TLAI', 'SOIL3C', 'TOTSOMC', 'TOTSOMC_1m', 'LEAFC', \
                       'DEADSTEMC', 'DEADCROOTC', 'FROOTC', 'LIVESTEMC', 'LIVECROOTC', 'TOTVEGC', \
                       'TOTCOLC', 'TOTLITC', 'BTRAN', 'SCALARAVG_vr', 'CWDC', 'QVEGE', 'QVEGT', 'QSOIL', 'QDRAI', \
                       'QRUNOFF', 'FPI', 'FPI_vr', 'FPG', 'FPI_P','FPI_P_vr', 'FPG_P', 'CPOOL','NPOOL', 'PPOOL', 'SMINN', 'HR_vr']
    if (options.C14):
        var_list_spinup.append('C14_TOTSOMC')
        var_list_spinup.append('C14_TOTSOMC_1m')
        var_list_spinup.append('C14_TOTVEGC')
    #ILAMB diagnostic variables
    ilamb_outputs = ['FAREA_BURNED', 'CWDC', 'LEAFC', 'TOTLITC', 'STORVEGC', 'LIVESTEMC', 'DEADSTEMC', \
                     'TOTPRODC', 'FROOTC', 'LIVECROOTC', 'DEADCROOTC', 'SOIL1C', 'SOIL2C', 'SOIL3C', \
                     'TOTSOMC', 'TOTVEGC', 'WOODC', 'QSOIL', 'QVEGE', 'COL_FIRE_CLOSS', \
                     'LITR1C_TO_SOIL1C', 'LITR2C_TO_SOIL2C', 'LITR3C_TO_SOIL3C', 'LAND_USE_FLUX', \
                     'LITFALL', 'GPP', 'FGR', 'TLAI', 'SNOWLIQ', 'SOILICE', 'SOILLIQ', 'QRUNOFF', \
                     'QOVER', 'SOILWATER_10CM', 'NBP', 'LEAFC_ALLOC', 'WOODC_ALLOC', 'QINTR', \
                     'AR', 'GR', 'HR', 'MR', 'FSNO', 'SNOWDP', 'QMELT', 'H2OSNO', 'SNOWBCMSL', \
                     'SNODSTMSL', 'SNOOCMSL', 'QVEGT', 'TSOI', 'WIND', 'EFLX_LH_TOT', 'FCTR', \
                     'FCEV', 'FGEV', 'FSH', 'RH2M', 'Q2M', 'RAIN', 'SNOW', 'PBOT', 'FLDS', 'FIRE', \
                     'FSDS', 'FSR', 'TSA']
    if ('CTC' in compset):
        var_list_daily.append('SOIL4C_vr')
        var_list_spinup.append('SOIL4C')
        ilamb_outputs.append('SOIL4C')

    if ('20TR' not in compset and int(options.hist_mfilt) == -1):
	#default to annual for spinup runs if not specified
	options.hist_mfilt = 1
	options.hist_nhtfrq = -8760

    if (options.hist_mfilt != -1 and not options.diags):
        if (options.ad_spinup):
            output.write(" hist_mfilt = "+str(options.hist_mfilt)+", "+str(options.hist_mfilt)+"\n")
        else:
            output.write(" hist_mfilt = "+ str(options.hist_mfilt)+"\n")
    if (options.hist_nhtfrq != -999 and not options.diags):
        if (options.ad_spinup):
            output.write(" hist_nhtfrq = "+ str(options.hist_nhtfrq)+", "+str(options.hist_nhtfrq)+"\n")
        else:
            output.write(" hist_nhtfrq = "+ str(options.hist_nhtfrq)+"\n")
    
    if (options.hist_vars != ''):
        output.write(" hist_empty_htapes = .true.\n")
        #read hist_vars file
        hvars_file = open(hist_vars)
        myline = " hist_fincl1 = "
        line2 = 0
        for s2 in hvars_file:
            if line2 ==0:
                myline = myline+"'"+s2.strip()+"'"
            else:
                myline = myline+",'"+s2.strip()+"'"
            line2=line2+1
        output.write(myline+"\n")
        hvars_file.close()
   
    if (options.spinup_vars and (not '20TR' in compset)):
        output.write(" hist_empty_htapes = .true.\n")
        h0varst = " hist_fincl1 = "
        for v in var_list_spinup:
	    h0varst = h0varst+"'"+v+"',"
        h0varst = h0varst[:-1]+"\n"
        output.write(h0varst)

    if ('20TR' in compset and options.diags):
        output.write(" hist_dov2xy = .true., .true., .true., .false., .true.\n")
        output.write(" hist_empty_htapes = .true.\n")
        output.write(" hist_mfilt = 1, 8760, 365, 365, 1\n")
        output.write(" hist_nhtfrq = 0, -1, -24, -24, -8760\n")
        h1varst = " hist_fincl2 = "
        h2varst = " hist_fincl3 = "
        h3varst = " hist_fincl4 = "
        h4varst = " hist_fincl5 = "
        for v in var_list_hourly:
            h1varst = h1varst+"'"+v+"',"
            h2varst = h2varst+"'"+v+"',"
            h4varst = h4varst+"'"+v+"',"
        for v in var_list_daily:
            h2varst = h2varst+"'"+v+"',"
            h4varst = h4varst+"'"+v+"',"
        for v in var_list_pft:
            h3varst = h3varst+"'"+v+"',"
        h1varst = h1varst[:-1]+"\n"
        h2varst = h2varst[:-1]+"\n"
        h3varst = h3varst[:-1]+"\n"
        h4varst = h4varst[:-1]+"\n"
        output.write(h1varst)
        output.write(h2varst)
        output.write(h3varst)
        output.write(h4varst)
    elif ('20TR' in compset and (options.trans_varlist != '' or options.ilambvars)):
	trans_varlist = options.trans_varlist.split(',')
        if (options.ilambvars):
            trans_varlist = ilamb_outputs
	output.write(" hist_empty_htapes = .true.\n")
        h0varst = " hist_fincl1 = "
	for v in trans_varlist:
	    h0varst = h0varst+"'"+v+"',"
	h0varst = h0varst[:-1]+"\n"
        output.write(h0varst)

    if (options.ad_spinup):
        output.write(" hist_dov2xy = .true., .false.\n")
        if ('CNT' in compset):
            output.write(" hist_fincl2 = 'CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL3C_vr', 'SOIL3N_vr', " \
                             +"'SOIL3P_vr', 'SOIL2C_vr','SOIL2N_vr', 'SOIL2P_vr', 'DEADSTEMC', " \
                             +"'DEADSTEMN', 'DEADSTEMP', 'DEADCROOTC', 'DEADCROOTN', "+ \
                             "'DEADCROOTP', 'LITR3C_vr', 'LITR3N_vr', 'LITR3P_vr'\n")
        else:	
            output.write(" hist_fincl2 = 'CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL4C_vr', 'SOIL4N_vr', " \
                             +"'SOIL4P_vr', 'SOIL3C_vr', 'SOIL3N_vr', 'SOIL3P_vr', 'DEADSTEMC', " \
                             +"'DEADSTEMN', 'DEADSTEMP', 'DEADCROOTC', 'DEADCROOTN', "+ \
                             "'DEADCROOTP', 'LITR3C_vr', 'LITR3N_vr', 'LITR3P_vr'\n")
        output.write(" finidat = ''\n")
    elif (options.coldstart == False):
        #user-defined initial data file
        output.write(" finidat = '"+finidat+"'\n")
        
    #surface data file
    output.write(" fsurdat = '"+rundir+"/surfdata.nc'\n")
        
    #pft dynamics file for transient run
    if ('20TR' in compset):
        if (options.nopftdyn):
            output.write(" flanduse_timeseries = ' '\n") 
        else:
            output.write(" flanduse_timeseries = '"+rundir+"/surfdata.pftdyn.nc'\n")
        output.write(' check_finidat_fsurdat_consistency = .false.\n')
        output.write(' check_finidat_year_consistency = .false.\n')
    #pft-physiology file
    output.write(" paramfile = '"+rundir+"/clm_params.nc'\n")


    if ('RD' in compset or 'ECA' in compset):
        #soil order parameter file
        output.write(" fsoilordercon = '"+rundir+"/CNP_parameters.nc'\n")
        output.write( " stream_fldfilename_ndep = '"+options.ccsm_input+ \
          "/lnd/clm2/ndepdata/fndep_clm_hist_simyr1849-2006_1.9x2.5_" + \
                      "c100428.nc'\n")
        if (options.vsoilc):
            output.write(" use_vertsoilc = .true.\n")
        if (options.centbgc):
            output.write(" use_century_decomp = .true.\n")
        if (options.no_dynroot):
            output.write(" use_dynroot = .false.\n")
        if (options.bulk_denitrif):
            output.write(" use_nitrif_denitrif = .false.\n")
        else:
            output.write(" use_nitrif_denitrif = .true.\n")
        if (options.CH4 or (not options.bulk_denitrif)):
            output.write(" use_lch4 = .true.\n")
        if (options.nofire):
            output.write(" use_nofire = .true.\n")
        if (options.c_only):
            options.write(" suplphos = 'ALL'\n")
            options.write(" suplnitro  = 'ALL'\n")
        elif (options.cn_only or options.ad_spinup):
            output.write(" suplphos = 'ALL'\n")
            output.write(" suplnitro = 'NONE'\n")
        elif (options.cp_only):
            output.write(" suplphos = 'NONE'\n")
            output.write(" suplnitro = 'ALL'\n")
        else:
            output.write(" suplphos = 'NONE'\n")
            output.write(" suplnitro = 'NONE'\n")
        if (options.C13):
            output.write(" use_c13 = .true.\n")
        if (options.C14):
            output.write(" use_c14 = .true.\n")
            output.write(" use_c14_bombspike = .true.\n")
            output.write(" atm_c14_filename = '"+options.ccsm_input+"/atm/datm7/CO2/" + \
                         "atm_delta_C14_data_1850-2007_monthly_25082011.nc'\n")
        if ('ECA' in compset):
            output.write(" nyears_ad_carbon_only = 0\n")
            output.write(" spinup_mortality_factor = 1\n")
        else:
            output.write(" nyears_ad_carbon_only = 25\n")
            output.write(" spinup_mortality_factor = 10\n")
    if (cpl_bypass):
        if (use_reanalysis):
            if (options.livneh):
                output.write(" metdata_type = 'livneh'\n")
                output.write(" metdata_bypass = '"+options.ccsm_input+"/atm/datm7/" \
                         +"atm_forcing.datm7.cruncep_qianFill.0.5d.V5.c140715_Livneh" + \
                         "/cpl_bypass_full'\n")
            elif (options.cruncep):
                output.write(" metdata_type = 'cru-ncep'\n")
                output.write(" metdata_bypass = '"+options.ccsm_input+"/atm/datm7/" \
                         +"atm_forcing.datm7.cruncep_qianFill.0.5d.V5.c140715/cpl_bypass_full'\n")
#                         +"atm_forcing.datm7.cruncep.0.5d._v4_c110920.ornl/cpl_bypass_full'\n")
            elif (options.gswp3):
                output.write(" metdata_type = 'gswp3'\n")
                output.write(" metdata_bypass = '"+options.ccsm_input+"/atm/datm7/" \
                         +"atm_forcing.datm7.GSWP3.0.5d.v1.c170516/cpl_bypass_full'\n")
            elif (options.cplhist):
                output.write(" metdata_type = 'cplhist'\n")
                output.write(" metdata_bypass = '"+options.ccsm_input+"/atm/datm7/" \
                         +"atm_forcing.cpl.WCYCL1850S.ne30.c171204/cpl_bypass_full'\n")
        else:
            output.write("metdata_type = 'site'\n")
            output.write(" metdata_bypass = '"+options.ccsm_input+"/atm/datm7/" \
                             +"CLM1PT_data/"+ptstr+"_"+options.site+"/'\n")
        output.write(" co2_file = '"+options.ccsm_input+"/atm/datm7/CO2/" \
                         +options.co2_file+"'\n")
        output.write(" aero_file = '"+options.ccsm_input+"/atm/cam/chem/" \
                         +"trop_mozart_aero/aero/aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'\n")


    output.close()

#configure case
if (options.no_config == False):
    os.system('./case.setup')
else:
    print("Warning:  No case configure performed")
    sys.exit()

#stream file modificaitons: directory and domain file (for using site_level CRU-NCEP)
if (not cpl_bypass and not isglobal):
    if (options.cruncep):
        types = ['Precip', 'Solar', 'TPQW']
        tout  = ['Precip', 'Solar', 'TPHWL']
        for i in range(0,3):
            input  = open('./Buildconf/datmconf/datm.streams.txt.CLMCRUNCEP.'+types[i], 'r')
            output = open('./user_datm.streams.txt.CLMCRUNCEP.'+types[i],'w')
            line = 1
            for s in input:
                if (line == 16):
                    output.write('            domain.lnd.'+ptstr+'_'+casename+'_navy.nc\n')
                elif (i < 2 and line == 24):
                    output.write('            '+options.ccsm_input+'/atm/datm7/ugrid/'+ptstr+'_'+options.site+ \
                                     '/'+tout[i]+'6Hrly\n')
                elif (i == 2 and line == 27):
                    output.write('            '+options.ccsm_input+'/atm/datm7/ugrid/'+ptstr+'_'+options.site+ \
                                     '/'+tout[i]+'6Hrly\n')
                else:
                    output.write(s)
                line = line+1
            input.close()
            output.close()
    if ('1850' in compset):
        myinput  = open('./Buildconf/datmconf/datm.streams.txt.presaero.clim_1850')
        myoutput = open('./user_datm.streams.txt.presaero.clim_1850','w')
        for s in myinput:
            if ('aerosoldep_monthly' in s):
                myoutput.write('            aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc\n')
            else:
                myoutput.write(s)
        myinput.close()
        myoutput.close()

    #reverse directories for CLM1PT and site
    if (options.cruncep == False):
        myinput  = open('./Buildconf/datmconf/datm.streams.txt.CLM1PT.CLM_USRDAT')
        myoutput = open('./user_datm.streams.txt.CLM1PT.CLM_USRDAT','w')
        for s in myinput:
            if ('CLM1PT_data' in s):
                temp = s.replace('CLM1PT_data', 'TEMPSTRING')
                s    = temp.replace(str(numxpts)+'x'+str(numypts)+'pt'+'_'+options.site, 'CLM1PT_data')
                temp  =s.replace('TEMPSTRING', str(numxpts)+'x'+str(numypts)+'pt'+'_'+options.site)
                myoutput.write(temp)
            else:
                myoutput.write(s)
        myinput.close()
        myoutput.close()

#CPPDEF modifications
infile  = open("./Macros.make")
outfile = open("./Macros.make.tmp",'a')
for s in infile:
    if (s[0:7] == "CPPDEFS"):
        stemp = s
        if (options.harvmod):
            print("Turning on HARVMOD modificaiton\n")
            stemp = stemp[:-1]+' -DHARVMOD\n'
        if (cpl_bypass):
            stemp = stemp[:-1]+' -DCPL_BYPASS\n'
        if (options.humhol):
            stemp = stemp[:-1]+' -DHUM_HOL\n'
        if (options.MARSH):
            stemp = stemp[:-1]+' -DMARSH\n'
        outfile.write(stemp) 
    elif (s[0:13] == "NETCDF_PATH:=" and options.machine == 'userdefined'):
        try:
            os.environ['NETCDF_PATH']
        except KeyError:
            print('ERROR:  Must set NETCDF_PATH environment variable for user defined machine')
            sys.exit(1)
        outfile.write('NETCDF_PATH:= '+os.getenv('NETCDF_PATH')+'\n')
    elif (s[0:7] == 'SLIBS+=' and options.machine == 'userdefined'):
        outfile.write('SLIBS+=-lnetcdff -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas' \
                          +' $(shell $(NETCDF_PATH)/bin/nc-config --flibs)\n')
    else:
        outfile.write(s)
infile.close()
outfile.close()
os.system('mv Macros.make.tmp Macros.make')


#copy sourcemods
os.chdir('..')
if (options.srcmods_loc != ''):
    if (os.path.exists(options.srcmods_loc) == False):
        print('Invalid srcmods directory.  Exiting')
        sys.exit()
    options.srcmods_loc = os.path.abspath(options.srcmods_loc)
    os.system('cp -r '+options.srcmods_loc+'/* ./'+casename+ \
                  '/SourceMods')
if (options.caseroot == './' ):
    os.chdir(csmdir+"/cime/scripts/"+casename)
else:
    os.chdir(casedir)       

#Datm mods/ transient CO2 patch for transient run (datm buildnml mods)
if (not cpl_bypass):
    myinput  = open('./Buildconf/datmconf/datm_in')
    myoutput = open('user_nl_datm','w')
    for s in myinput:
        if ('streams =' in s):
            myalign_year = 1 #startyear
            if (options.align_year != -999):
                myalign_year = options.align_year
            if ('20TR' in compset):
                mypresaero = '"datm.streams.txt.presaero.trans_1850-2000 1850 1850 2000"'
                myco2      = ', "datm.streams.txt.co2tseries.20tr 1766 1766 2010"'
            else:
                mypresaero = '"datm.streams.txt.presaero.clim_1850 1 1850 1850"'
                myco2=''
            if (options.cruncep):
                myoutput.write(' streams = "datm.streams.txt.CLMCRUNCEP.Solar '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+ \
                                   '"datm.streams.txt.CLMCRUNCEP.Precip '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+ \
                                   '"datm.streams.txt.CLMCRUNCEP.TPQW '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+mypresaero+myco2+ \
                                   ', "datm.streams.txt.topo.observed 1 1 1"\n')
            else:
                myoutput.write(' streams = "datm.streams.txt.CLM1PT.CLM_USRDAT '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+mypresaero+myco2+ \
                                   ', "datm.streams.txt.topo.observed 1 1 1"\n')
        elif ('streams' in s):
            continue  #do nothing
        elif ('taxmode' in s):
            if (options.cruncep):
                taxst = "taxmode = 'cycle', 'cycle', 'cycle', 'extend', 'extend'"
            else:
                taxst = "taxmode = 'cycle', 'extend', 'extend'"
            if ('20TR' in compset):
                taxst = taxst+", 'extend'"
            myoutput.write(taxst+'\n')
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()

os.system('mkdir Srcfiles')
#clean build if requested
if (options.clean_build):
    os.system('./case.build --clean')
#compile cesm
if (options.no_build == False):
    os.system('./case.build')
    if ('eos' in options.machine):
      os.system('./case.build')
else:
    print ('no_build set.  Assuming model has already been built.')
    print ('Creating run directory.  Setting BUILD_COMPLETE = TRUE')
    os.system('./xmlchange BUILD_COMPLETE=TRUE')
    os.system('mkdir -p '+rundir)
if (options.caseroot == ''):
    os.chdir(csmdir+"/cime/scripts/"+casename)
else:
    os.chdir(casedir) 

#copy site data to run directory
os.system('cp '+PTCLMdir+'/temp/domain.nc '+PTCLMdir+'/temp/surfdata.nc  '+ \
              PTCLMdir+'/temp/*param*.nc '+runroot+'/'+casename+'/run/')
if ('20TR' in compset and options.nopftdyn == False):
    os.system('cp '+PTCLMdir+'/temp/surfdata.pftdyn.nc '+runroot+'/'+casename+'/run/')

#submit job if requested
if (options.no_submit == False and int(options.mc_ensemble) < 0 and options.ensemble_file == ''):
# Ming, 01/22/16, use a csh script instead of the perl script to submit the job.
#    os.system("qsub "+casename+".run")
    if ("mesabi" in options.machine):
        os.system("cp /home/reichpb/shared/Scripts/Ming/shell/qsub1.csh ./")
        os.system("cp /home/reichpb/shared/Scripts/Ming/shell/st_archive.sh ./")
        f = open("qsub.csh",'r+')
        filedata = f.read()
        newdata0 = filedata.replace("newdir",casedir)
        newdata = newdata0.replace("archivedir",options.archiveroot)
        f.seek(0)
        f.write(newdata)
        f.truncate()
        f.close()
        os.system("qsub qsub.csh")
    else:
        os.system("./case.submit")



#------------------------- Code to generate and run parameter ensembles --------------------------------------

os.chdir(PTCLMdir)

if (options.ensemble_file != '' or int(options.mc_ensemble) != -1):
    if (not(os.path.isfile(options.parm_list))):
	print('parm_list file does not exist')
        sys.exit()
    elif (isglobal):
        print('Ensemble simulations not supported for regional/global simulations')
        sys.exit()
    else:
        param_names=[]
        param_min=[]
        param_max=[]
        input = open(options.parm_list,'r')
        for s in input:
	    if (s):
                param_names.append(s.split()[0])
                if (int(options.mc_ensemble) > 0):
                    if (len(s.split()) == 3):
                        param_min.append(float(s.split()[1]))
                        param_max.append(float(s.split()[2]))
                    else:
                        param_min.append(float(s.split()[2]))
                        param_max.append(float(s.split()[3]))
        input.close() 
        n_parameters = len(param_names)

    if (options.ensemble_file != ''):    
        if (not os.path.isfile(options.ensemble_file)):
            print('Error:  ensemble file does not exist')
            sys.exit()

        samples=numpy.zeros((n_parameters,100000), dtype=numpy.float) 
        #get parameter samples and information
        myinput=open(options.ensemble_file)
        nsamples = 0
        for s in myinput:
            for j in range(0,n_parameters):
                samples[j][nsamples] = float(s.split()[j]) 
            nsamples=nsamples+1
        myinput.close()
    elif (int(options.mc_ensemble) > 0):
        nsamples = int(options.mc_ensemble)
        samples=numpy.zeros((n_parameters,nsamples), dtype=numpy.float)
        for i in range(0,nsamples):
            for j in range(0,n_parameters):
                samples[j][i] = param_min[j]+(param_max[j]-param_min[j])*numpy.random.rand(1)
        numpy.savetxt('mcsamples_'+casename+'.txt', numpy.transpose(samples))
        options.ensemble_file = 'mcsamples_'+casename+'.txt'

    print str(n_parameters)+' parameters are being modified' 
    print str(nsamples)+' parameter samples provided'
  
    #total number of processors required in each pbs script
    np_total = int(options.np)*int(options.ng)
    #number of scripts required
    n_scripts = int(math.ceil(nsamples/float(options.ninst*options.ng)))
 
    num=0
    #Launch ensemble if requested 
    mysubmit_type = 'qsub'
    if ('cori' in options.machine or options.machine == 'edison'):
        mysubmit_type = 'sbatch'
    if (options.ensemble_file != ''):
        #write the PBS scripts to copy directories and run the ensembles for this case
        if (os.path.exsits(casedir+'/case.run')):
            myinput = open(casedir+'/case.run','r')
        else:
            myinput = open(casedir+'./case.run','r')
        output_run  = open(tmpdir+'/ensemble_run_'+casename+'.pbs','w')
        timestr=str(int(float(options.walltime)))+':'+str(int((float(options.walltime)- \
                                     int(float(options.walltime)))*60))+':00'
        for s in myinput:
            if ("perl" in s or 'python' in s):
                output_run.write("#!/bin/csh -f\n")
                if (mysubmit_type == 'qsub'):
                    output_run.write('#PBS -l walltime='+timestr+'\n')
                else:
                    output_run.write('#SBATCH --time='+timestr+'\n')
                    if ('edison' in options.machine or 'cori' in options.machine):
                        if (options.debug):
                            output_run.write('#SBATCH --partition=debug\n')
                        else:
	                    output_run.write('#SBATCH --partition=regular\n')
            elif ("#PBS" in s or "#!" in s or '#SBATCH' in s):
                #edit number of required nodes for ensemble runs
                if ('nodes' in s):
                    if (options.machine == 'cades'):
                        output_run.write('#PBS -l nodes='+str(int(math.ceil(np_total/(ppn*1.0))))+ \
                                             ':ppn='+str(ppn)+'\n')
                    elif (mysubmit_type == 'qsub'):
                        output_run.write('#PBS -l nodes='+str(int(math.ceil(np_total/(ppn*1.0))))+ \
                                         '\n')
                    else:
                        output_run.write('#SBATCH --nodes='+str(int(math.ceil(np_total/(ppn*1.0))))+ \
                                         '\n')
                else:
                    output_run.write(s)
        myinput.close()

        output_run.write("\n")
        if (options.machine == 'eos'):
            output_run.write('source $MODULESHOME/init/csh\n')
            output_run.write('module load nco\n')
            output_run.write('module load cray-netcdf\n')
	    output_run.write('module unload python\n')
            output_run.write('module load python/2.7.5\n')
            output_run.write('module unload PrgEnv-intel\n')
            output_run.write('module load PrgEnv-gnu\n')
            output_run.write('module load python_numpy\n')
            output_run.write('module load python_scipy\n')
            output_run.write('module load python_mpi4py/2.0.0\n')
            output_run.write('module unload PrgEnv-gnu\n')
            output_run.write('module load PrgEnv-intel\n')
        if (options.machine == 'titan'):
            output_run.write('source $MODULESHOME/init/csh\n')
            output_run.write('module load nco\n')
            output_run.write('module load cray-netcdf\n')
            output_run.write('module load python/2.7.9\n')
            output_run.write('module load python_numpy/1.9.2\n')
            output_run.write('module load python_scipy/0.15.1\n')
            output_run.write('module load python_mpi4py/2.0.0\n')
        output_run.write('cd '+csmdir+'/components/clm/tools/clm4_5/pointclm/\n')
        cnp = 'True'
        if (options.cn_only or options.c_only):
            cnp= 'False'
        if ('oic' in options.machine or 'cades' in options.machine):
            mpicmd = 'mpirun'
            if ('cades' in options.machine):
                mpicmd = '/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/1.10.3/centos7.2_gnu5.3.0/bin/mpirun'
            cmd = mpicmd+' -np '+str(np_total)+' --hostfile $PBS_NODEFILE python manage_ensemble.py ' \
               +'--case '+casename+' --runroot '+runroot+' --n_ensemble '+str(nsamples)+' --ens_file '+ \
               options.ensemble_file+' --exeroot '+exeroot+' --parm_list '+options.parm_list+' --cnp '+cnp + \
               ' --site '+options.site
        elif (('titan' in options.machine or 'eos' in options.machine or 'edison' in options.machine) and int(options.ninst) == 1):
            cmd = 'aprun -n '+str(np_total)+' python manage_ensemble.py ' \
               +'--case '+casename+' --runroot '+runroot+' --n_ensemble '+str(nsamples)+' --ens_file '+ \
               options.ensemble_file+' --exeroot '+exeroot+' --parm_list '+options.parm_list+' --cnp '+cnp + \
               ' --site '+options.site
        elif ('cori' in options.machine):
            cmd = 'srun -n '+str(np_total)+' python manage_ensemble.py ' \
               +'--case '+casename+' --runroot '+runroot+' --n_ensemble '+str(nsamples)+' --ens_file '+ \
               options.ensemble_file+' --exeroot '+exeroot+' --parm_list '+options.parm_list+' --cnp '+cnp + \
               ' --site '+options.site

        if (options.postproc_file != ''): 
            cmd = cmd + ' --postproc_file '+options.postproc_file
        output_run.write(cmd+'\n')
        output_run.close()
        if (options.no_submit == False):
            os.system('qsub '+tmpdir+'/ensemble_run_'+casename+'.pbs')
