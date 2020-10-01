#!/usr/bin/env python

import netcdf_functions as nffun
import os, sys, csv, time, math, numpy
import subprocess
from optparse import OptionParser
#from Numeric import *


#DMR 4/16/13
#call_runCLM.py does the following:
#  1. Call routines to create point data (makepointdata.py, makemetdata.py)
#  2. Set point and case-specific namelist options
#  2. configure case
#  3. build (compile) CESM with clean_build first if requested
#  4. apply patch for transient CO2 if transient run
#  6. apply user-specified PBS and submit information
#  7. submit job to PBS queue if requested.
#
#  For reproducibility, a copy of the current call_PTCLM.py is saved
#  to the newly created case directory.  This is for informational
#  purposes only - the script should not be executed from within
#  the case directory.
#
# Add changes from FMY 6/6/2013
# modified to work for CLM4-pf (CLM4.5.10, with PFLOTRAN interface) version used by NGEE-Arc



#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--runroot", dest="runroot", default="", \
                  help="Directory where the run would be created")
parser.add_option("--exeroot", dest="exeroot", default="", \
	         help="Location of executable")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--coldstart", dest="coldstart", default=False, \
                  help = "set cold start (mutually exclusive w/finidat)", \
                  action="store_true")
parser.add_option("--compset", dest="compset", default='I1850CLM45CN', \
                  help = "component set to use (required)\n"
                         "Currently supports ONLY *CLM45(CN) compsets")
parser.add_option("--cruncep", dest="cruncep", default=False, \
                  help = "use cru-ncep data", action="store_true")
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
parser.add_option("--regional", action="store_true", \
                   dest="regional", default=False, \
                   help="Flag for regional run (2x2 or greater)")
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
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=600, \
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
parser.add_option("--xpts", dest="xpts", default=1, \
                      help = 'for regional runs: xpts')
parser.add_option("--ypts", dest="ypts", default=1, \
                      help = 'for regional runs: ypts')
parser.add_option("--trans2", dest="trans2", default=False, action="store_true", \
                  help = 'Tranisnent phase 2 (1901-2010) - CRUNCEP only')
parser.add_option("--spinup_vars", dest="spinup_vars", default=False, \
                  help = 'Limit output vars in spinup runs', action="store_true")
parser.add_option("--cn_only", dest="cn_only", default=True, \
                  help = 'Carbon/Nitrogen only (supplemental P)', action="store_true")
parser.add_option("--cnp", dest="cnp", default=False, \
                  help = 'CNP model', action = "store_true")
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

(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------

#machine info:  cores per node
ppn=1
if ('titan' in options.machine):
    ppn=16
elif ('oic2' in options.machine):
    ppn=8
elif ('oic5' or 'eos' or 'cades' in options.machine):
    ppn=32
elif ('edison' in options.machine):
    ppn=24

print ppn
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

if ('CB' in compset):
    cpl_bypass = True
else:
    cpl_bypass = False

#figure out if clm40 or clm45, set model-specific options
isclm45 =  False
surfdir = 'surfdata'
pftphys_stamp = 'clm40.c130424'
if ('CLM45' in compset):
    isclm45 = True
    surfdir = 'surfdata_map'
    pftphys_stamp = 'c160128'
CNPstamp = 'c131108'


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

if ('CN' in compset):
  mybgc = 'CN'
elif ('BGC' in compset):
  mybgc = 'BGC'
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
if ('CRU' in compset or options.cruncep):
    use_cruncep = True
else:
    use_cruncep = False
if ('20TR' in compset and use_cruncep and not ('CB' in compset)):
    if options.trans2:
        casename = casename+'_phase2'
    else:
        casename = casename+'_phase1'
if (options.ad_spinup):
    casename = casename+'_ad_spinup'
if (options.exit_spinup):
    casename = casename+'_exit_spinup'

#PTCLMfiledir = csmdir+'/scripts/acme/pointclm/PTCLM_files/'
PTCLMfiledir = options.ccsm_input+'/lnd/clm2/PTCLM'

if (caseroot != "./"):
    casedir=caseroot+"/"+casename
else:
    casedir=casename

#Check for existing case directory
if (os.path.exists(casedir)):
    
    print('Warning:  Case directory exists')
    if (options.rmold):
        print('--rmold specified.  Removing old case (this will NOT clean the run directory')
        print('Please perform a clean build if code has changed')
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
    if ('titan' in options.machine):
        exeroot = os.path.abspath(os.environ['HOME']+ \
    	    '/acme_scratch/pointclm/'+casename+'/bld')
else:
    exeroot=options.exeroot
print("CASE exeroot is: "+exeroot+"\n")
rundir=runroot+'/'+casename+'/run'
print("CASE rundir is: "+rundir+"\n")


#------------------- make point data for site -------------------------------
if (options.nopointdata == False):
    ptcmd = 'python makepointdata.py --caseroot '+caseroot+' --casename '+casename+ \
        ' --site '+options.site+' --sitegroup '+options.sitegroup+ \
        ' --csmdir '+csmdir+' --ccsm_input '+options.ccsm_input+ \
        ' --compset '+compset

    if (options.metdir != 'none'):
        ptcmd = ptcmd + ' --metdir '+options.metdir
    if (options.makemet):
        ptcmd = ptcmd + ' --makemetdata'
    if (options.surfdata_grid):
        ptcmd = ptcmd + ' --surfdata_grid'
    if (options.include_nonveg):
        ptcmd = ptcmd + ' --include_nonveg'
    if (options.regional):
        ptcmd = ptcmd + ' --regional'
        ptcmd = ptcmd + ' --xpts '+options.xpts
        ptcmd = ptcmd + ' --ypts '+options.ypts
    if ('45' not in compset):
        ptcmd = ptcmd + ' --clm40'
    print(ptcmd)
    os.system(ptcmd)
else:
    print('point data making NOT requested!  Make sure they exist')

#get site year information
sitedatadir = os.path.abspath(PTCLMfiledir)
os.chdir(sitedatadir)
AFdatareader = csv.reader(open(options.sitegroup+'_sitedata.txt',"rb"))
for row in AFdatareader:
    if row[0] == options.site:
        if (use_cruncep):
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
        if (options.regional == True):
            if (options.xpts < 2 and options.ypts < 2):
                print('Error:  xpts OR ypts MUST be greater than 1 for regional option\n')
                sys.exit()
            numxpts = int(options.xpts)
            numypts = int(options.ypts)
        else:
            numxpts=1
            numypts=1

           #numxpts=int(row[9])
           #numypts=int(row[10])
ptstr = str(numxpts)+'x'+str(numypts)+'pt'
os.chdir(csmdir+'/cime/scripts')
#get simyr
mysimyr=1850
if (options.compset == 'ICLM45CN' or options.compset == 'ICLM45BGC' or '2000' in compset):
    mysimyr=2000

#parameter (pft-phys) modifications if desired
tmpdir = csmdir+'/components/clm/tools/clm4_5/pointclm/temp'
os.system('mkdir -p '+tmpdir)
os.system('cp '+options.ccsm_input+'/lnd/clm2/paramdata/clm_params.'+pftphys_stamp+'.nc ' \
              +tmpdir+'/clm_params.'+pftphys_stamp+'.'+casename+'.nc')
os.system('chmod u+w ' +tmpdir+'/clm_params.'+pftphys_stamp+'.'+casename+'.nc')
if (options.parm_file != ''):
    pftfile = tmpdir+'/clm_params.'+pftphys_stamp+'.'+casename+'.nc'
    input   = open(os.path.abspath(options.parm_file))
    for s in input:
        if s[0:1] != '#':
            values = s.split()
            thisvar = nffun.getvar(pftfile, values[0])
            thisvar[int(values[1])] = float(values[2])
            ierr = nffun.putvar(pftfile, values[0], thisvar)
    input.close()

#parameter (soil order dependent) modifications if desired    ::X.YANG 
os.system('cp '+options.ccsm_input+'/lnd/clm2/paramdata/CNP_parameters_'+CNPstamp+'.nc ' \
              +tmpdir+'/CNP_parameters_'+CNPstamp+casename+'.nc')
os.system('chmod u+w ' +tmpdir+'/CNP_parameters_'+CNPstamp+casename+'.nc')
if (options.parm_file_P != ''):
    soilorderfile = options.ccsm_input+'/lnd/clm2/paramdata/CNP_parameters_'+CNPstamp+casename+'.nc'
    input   = open(os.path.abspath(options.parm_file_P))
    for s in input:
        if s[0:1] != '#':
            values = s.split()
            thisvar = nffun.getvar(pftfile, values[0])
            thisvar[int(values[1])] = float(values[2])
            ierr = nffun.putvar(pftfile, values[0], thisvar)
    input.close()
    soilorderfile.close()

#set number of run years for ad, exit spinup cases
if (options.ny_ad != options.run_n and options.ad_spinup):
    options.run_n = options.ny_ad
elif (options.exit_spinup):
    options.run_n = 1


#create new case
print ('./create_newcase --case '+casename+' --mach '+options.machine+' --compset '+ \
           options.compset+' --res CLM_USRDAT --compiler '+options.compiler+' --mpilib '+ \
           options.mpilib)
os.system('./create_newcase --case '+casename+' --mach '+options.machine+' --compset '+ \
	   options.compset+' --res CLM_USRDAT --compiler '+options.compiler+' --mpilib '+ \
           options.mpilib+' > create_newcase.log')
if (os.path.isdir(casename)):
    print(casename+' created.  See create_newcase.log for details')
    os.system('mv create_newcase.log '+casename)
else:
    print('failed to create case.  See create_newcase.log for details')

os.chdir(casedir)

#------------------ env_build.xml modifications ------------------------
os.system('./xmlchange -file env_build.xml -id EXEROOT -val '+exeroot)

#turn off ROF module
os.system('./xmlchange -file env_build.xml -id RTM_MODE -val NULL')

#clm 4_5 cn config options
if (isclm45):
    clmcn_opts = "-phys clm4_5"
else:
    clmcn_opts = "-phys clm4_0"

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
os.system('./xmlchange -file env_run.xml -id CLM_USRDAT_NAME ' \
              +' -val '+str(numxpts)+'x'+str(numypts)+'pt_'+options.site)
if (options.ad_spinup):
    os.system('./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS ' \
                  +' -val "-mask navy -bgc_spinup on -bgc '+mybgc.lower()+'"')
elif ('CN' in compset or 'BGC' in compset):
    os.system('./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS ' \
                  +' -val "-mask navy -bgc '+mybgc.lower()+'"')
else:
    os.system('./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS ' \
                  +' -val "-mask navy"')
if (options.run_startyear > -1):
    os.system('./xmlchange -file env_run.xml -id RUN_STARTDATE ' \
                  +' -val '+str(options.run_startyear)+'-01-01')
os.system('./xmlchange -file env_run.xml -id ATM_DOMAIN_PATH ' \
                  +' -val "\${RUNDIR}"')
os.system('./xmlchange -file env_run.xml -id LND_DOMAIN_PATH ' \
                  +' -val "\${RUNDIR}"')
#turn off archiving
os.system('./xmlchange -file env_run.xml -id DOUT_S ' \
                  +' -val "FALSE"') 
#datm options
if (use_cruncep):
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_MODE -val CLMCRUNCEP') 
else:
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

#if number of land instances > 1
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val '+str(options.np))
os.system('./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE -val '+str(options.np))

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
    os.system('./cesm_setup -clean')
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
    if (options.hist_mfilt != -1):
        if (options.ad_spinup):
            output.write(" hist_mfilt = "+str(options.hist_mfilt)+", "+str(options.hist_mfilt)+"\n")
        else:
            output.write(" hist_mfilt = "+ str(options.hist_mfilt)+"\n")
    if (options.hist_nhtfrq != -999):
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
        print myline
        output.write(myline+"\n")
        hvars_file.close()
    if (options.spinup_vars and (not '20TR' in compset)):
        output.write(" hist_empty_htapes = .true.\n")
        output.write(" hist_fincl1 = 'NPOOL', 'RETRANSN', 'PCO2', 'PBOT', 'NDEP_TO_SMINN', 'OCDEP', 'BCDEP', 'COL_FIRE_CLOSS', 'HDM', 'LNFM', 'NEE', 'GPP', 'FPSN', 'AR', 'HR', 'MR', 'GR', 'ER', 'NPP', 'TLAI', 'TOTSOMC', 'LEAFC', 'DEADSTEMC', 'DEADCROOTC', 'FROOTC', 'LIVESTEMC', 'LIVECROOTC', 'TOTVEGC', 'TOTCOLC', 'TOTLITC', 'BTRAN', 'CWDC', 'QVEGE', 'QVEGT', 'QSOIL', 'QDRAI', 'QRUNOFF', 'FPI', 'FPG'\n")

    if (options.ad_spinup):
        output.write(" hist_dov2xy = .true., .false.\n")
        if ('BGC' in compset or options.centbgc):

            output.write(" hist_fincl2 = 'CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL3C_vr', 'SOIL3N_vr', 'SOIL3P_vr', 'SOIL2C_vr', " + \
                             "'SOIL2N_vr', 'SOIL2P_vr', 'DEADSTEMC', 'DEADSTEMN', 'DEADSTEMP', 'DEADCROOTC', 'DEADCROOTN', "+ \
                             "'DEADCROOTP', 'LITR3C_vr', 'LITR3N_vr', 'LITR3P_vr'\n")
        else:	
            output.write(" hist_fincl2 = 'CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL4C_vr', 'SOIL4N_vr', 'SOIL4P_vr', 'SOIL3C_vr', " + \
	                     "'SOIL3N_vr', 'SOIL3P_vr', 'DEADSTEMC', 'DEADSTEMN', 'DEADSTEMP', 'DEADCROOTC', 'DEADCROOTN', "+ \
                             "'DEADCROOTP', 'LITR3C_vr', 'LITR3N_vr', 'LITR3P_vr'\n")

    #user-defined initial data file
    if (finidat != ''):
        output.write(" finidat = '"+finidat+"'\n")
    #surface data file

    if (options.nopointdata):
        output.write(" fsurdat = '"+options.ccsm_input+"/lnd/clm2/"+surfdir+"/surfdata_"+str(numxpts)+'x'+ \
                         str(numypts)+"pt_"+options.site+"_simyr"+str(mysimyr)+".nc'\n")
    else:
        output.write(" fsurdat = './surfdata_"+str(numxpts)+'x'+ \
                         str(numypts)+"pt_"+casename+"_simyr"+str(mysimyr)+".nc'\n")
    #pft dynamics file for transient run
    if ('20TR' in compset):
        if (options.nopftdyn):
            output.write(" flanduse_timeseries = ' '\n") 
        else:
            if (options.nopointdata):
	        output.write(" flanduse_timeseries = '"+options.ccsm_input+"/lnd/clm2/"+surfdir+ \
                  "/surfdata.pftdyn_"+str(numxpts)+'x'+str(numypts)+"pt_"+options.site+".nc'\n")
            else:
                output.write(" flanduse_timeseries = './surfdata.pftdyn_"+str(numxpts)+'x' \
                             +str(numypts)+"pt_"+casename+".nc'\n")
        output.write(' check_finidat_fsurdat_consistency = .false.\n')
        output.write(' check_finidat_year_consistency = .false.\n')
    #pft-physiology file
    if ("mesabi" in options.machine):
        if (options.mod_parm_file != ''):
            output.write(" paramfile = '"+options.mod_parm_file+"'\n")
    else:
        output.write(" paramfile = './clm_params."+pftphys_stamp+"."+ \
                     casename+".nc'\n")
    #soil order parameter file
    output.write(" fsoilordercon = './CNP_parameters_"+CNPstamp+ \
                     casename+".nc'\n")

    #nitrogen deposition file
    if ('CN' in compset or 'BGC' in compset):
        output.write( " stream_fldfilename_ndep = '"+options.ccsm_input+ \
        "/lnd/clm2/ndepdata/fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc'\n")
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
    if (options.nofire and isclm45):
        output.write(" use_nofire = .true.\n")
    if (options.cn_only or options.ad_spinup):
        output.write(" suplphos = 'ALL'\n")
    elif (options.cnp):
        output.write(" suplphos = 'NONE'\n")
    if (options.C13):
        output.write(" use_c13 = .true.\n")
    if (options.C14):
        output.write(" use_c14 = .true.\n")
    if (cpl_bypass):
        if (use_cruncep):
            output.write(" metdata_type = 'cru-ncep'\n")
            output.write(" metdata_bypass = '"+options.ccsm_input+"/atm/datm7/" \
                             +"atm_forcing.datm7.cruncep.0.5d._v4_c110920.ornl/cpl_bypass_full'\n")
        else:
            output.write("metdata_type = 'site'\n")
            output.write(" metdata_bypass = '"+options.ccsm_input+"/atm/datm7/" \
                             +"CLM1PT_data/"+ptstr+"_"+options.site+"/'\n")
        output.write(" co2_file = '"+options.ccsm_input+"/atm/datm7/CO2/" \
                         +options.co2_file+"'\n")
        output.write(" aero_file = '"+options.ccsm_input+"/atm/cam/chem/" \
                         +"trop_mozart_aero/aero/aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'\n")
    output.write(" nyears_ad_carbon_only = 25\n")
    output.write(" spinup_mortality_factor = 10\n")

    output.close()

#configure case
if (options.no_config == False):
    os.system('./cesm_setup')
else:
    print("Warning:  No case configure performed")
    sys.exit()

#stream file modificaitons: directory and domain file (for using site_level CRU-NCEP)
if (not cpl_bypass):
    if (use_cruncep):
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
            if (s[0:22] == '            aerosoldep'):
                myoutput.write('            aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc\n')
            else:
                myoutput.write(s)
        myinput.close()
        myoutput.close()

    #reverse directories for CLM1PT and site
    if (use_cruncep == False):
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
infile  = open("Macros")
outfile = open("Macros.tmp",'a')
for s in infile:
    if (s[0:7] == "CPPDEFS"):
        stemp = s
        if (options.nofire and isclm45 == False):
            print("Turning off FIRE\n")
            stemp = stemp[:-1]+' -DNOFIRE\n'
        if (options.harvmod):
            print("Turning on HARVMOD modificaiton\n")
            stemp = stemp[:-1]+' -DHARVMOD\n'
        if (options.ad_spinup and isclm45 == False):
            print("Turning on AD_SPINUP (CLM4.0)")
            stemp = stemp[:-1]+' -DAD_SPINUP\n'
        if (options.exit_spinup and isclm45 == False): 
            print("Turning on EXIT_SPINUP (CLM4.0)")
            stemp = stemp[:-1]+' -DEXIT_SPINUP\n'
        if (cpl_bypass):
            stemp = stemp[:-1]+' -DCPL_BYPASS\n'
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
os.system('mv Macros.tmp Macros')


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
    myinput  = open('./Buildconf/datmconf/datm_atm_in')
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
            if (use_cruncep):
                myoutput.write(' streams = "datm.streams.txt.CLMCRUNCEP.Solar '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+ \
                                   '"datm.streams.txt.CLMCRUNCEP.Precip '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+ \
                                   '"datm.streams.txt.CLMCRUNCEP.TPQW '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+mypresaero+myco2+'\n')
            else:
                myoutput.write(' streams = "datm.streams.txt.CLM1PT.CLM_USRDAT '+str(myalign_year)+ \
                                   ' '+str(startyear)+' '+str(endyear)+'  ", '+mypresaero+myco2+'\n')
        elif ('streams' in s):
            continue  #do nothing
        elif ('taxmode' in s):
            if (use_cruncep):
                taxst = "taxmode = 'cycle', 'cycle', 'cycle', 'extend'"
            else:
                taxst = "taxmode = 'cycle', 'extend'"
            if ('20TR' in compset):
                taxst = taxst+", 'extend'"
            myoutput.write(taxst+'\n')
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()

#clean build if requested
if (options.clean_build):
    os.system('./'+casename+'.clean_build')
#compile cesm
if (options.no_build == False):
    os.system('./'+casename+'.build')
else:
    print ('no_build set.  Assuming model has already been built.')
    print ('Creating run directory.  Setting BUILD_COMPLETE = TRUE')
    os.system('./xmlchange BUILD_COMPLETE=TRUE')
    os.system('mkdir -p '+rundir)
if (options.caseroot == ''):
    os.chdir(csmdir+"/cime/scripts/"+casename)
else:
    os.chdir(casedir) 

#move site data to run directory
os.system('mv '+tmpdir+'/*'+casename+'*.nc '+runroot+'/'+casename+'/run/')
if (options.nopointdata == False):
   os.system('mv '+tmpdir+'/domain*'+options.site+'* '+runroot+'/'+casename+'/run/')
else:
   os.system('cp '+options.ccsm_input+'/share/domains/domain.clm/domain*'+options.site+'* ' \
	+runroot+'/'+casename+'/run/')
if ('titan' in options.machine):
   os.system('cp '+exeroot+'/cesm.exe '+csmdir+'/components/clm/tools/clm4_5/pointclm/temp/')
#os.system('cp -f ../microbepar_in ' +csmdir+'/run/'+casename+'/run/')

#submit job if requested
if (options.no_submit == False and options.mc_ensemble < 0 and options.ensemble_file == ''):
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
        os.system("qsub "+casename+".run")

    
#------------------------- Code to generate and run parameter ensembles --------------------------------------

os.chdir(PTCLMdir)

if (options.ensemble_file != '' or options.mc_ensemble != -1):
    if (not(os.path.isfile('parm_list'))):
	print('parm_list file does not exist')
        sys.exit()
    else:
        param_names=[]
        param_min=[]
        param_max=[]
        input = open('parm_list','r')
        for s in input:
	    if (s):
                param_names.append(s.split()[0])
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
    elif (options.mc_ensemble > 0):
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
    #create ensemble directories 
    if (options.ensemble_file != '' or options.mc_ensemble != -1):
        for i in range(0,nsamples):
            #write the PBS scripts to copy directories and run the ensembles for this case
            if (i % int(options.ng) == 0):
                input = open(caseroot+'/'+casename+'/'+casename+'.run')
                numst=str(1000+num)
                output_copy = open(tmpdir+'/ensemble_copy_'+casename+'_'+numst[1:]+'.pbs','w')
                output_run  = open(tmpdir+'/ensemble_run_'+casename+'_'+numst[1:]+'.pbs','w')
                for s in input:
                    if ("perl" in s):
                        output_copy.write("#!/bin/csh -f\n")
                        output_run.write("#!/bin/csh -f\n")
                    elif ("#PBS" in s or "#!" in s):
                        #edit number of required nodes for ensemble runs
                        if ('nodes' in s):
                            output_copy.write('#PBS -l nodes=1:ppn='+str(ppn)+'\n')
                            output_run.write('#PBS -l nodes='+str(int(math.ceil(np_total/(ppn*1.0))))+ \
                                             ':ppn='+str(ppn)+'\n')
                        elif ('walltime' in s): 
                            output_copy.write('#PBS -l walltime=3:00:00\n')
                            output_run.write('#PBS -l walltime=3:00:00\n')
                        else:
                            output_copy.write(s)
                            output_run.write(s)
                input.close()
                output_copy.write("\n")
                output_run.write("\n")
            ngst=str(100000+(i % int(options.ng))+1)
            if ('oic' in options.machine or 'cades' in options.machine):
                #need to distribute jobs to nodes manually on oic
                myline_end = int(options.np)*((i % int(options.ng))+1)*int(options.ninst)
                output_run.write('head -'+str(myline_end)+' $PBS_NODEFILE | tail -1 > '+tmpdir+ \
                              '/mynodefile'+ngst[1:]+'\n')
            est=str(100000+i+1)
            ens_dir  = runroot+'/UQ/'+casename+'/g'+est[1:]
            if (int(est)-100000 <= math.ceil(float(nsamples)/options.ninst)):
              output_copy.write('cd '+csmdir+'/components/clm/tools/clm4_5/pointclm/\n')
              output_copy.write('python ensemble_copy.py --case '+casename+' --runroot '+runroot \
                                    +' --ens_num '+str(i+1)+' --ens_file '+options.ensemble_file+' &\n')
              if ((i+1) % ppn == 0 or (i == nsamples-1)):
	        output_copy.write('wait\n')
	        
              output_copy.write('mkdir -p '+ens_dir+'/timing/checkpoints\n')
              if ('oic' in options.machine or 'cades' in options.machine):
                   mpi_args = ' -np '+str(options.np)+' --hostfile '+tmpdir+'/mynodefile'+ngst[1:]
                   if ('cades' in options.machine):        #NOTE - will work with np=1 only
                  	mpi_args = mpi_args+' --cpu-set '+str(i % ppn)
                   output_run.write('cd '+ens_dir+'\n')
                   output_run.write('mpirun '+mpi_args+' '+exeroot+'/cesm.exe > ccsm_log.txt &\n')
              elif (('titan' in options.machine or 'eos' in options.machine) and int(options.ninst) == 1):
                  #use wraprun utility on nccs to manage the ensemble
                  output_run.write('cd '+ens_dir+'\n')
                  if ( (i % ppn) == 0):
                      cmd = 'wraprun -n '
                      for pp in range(0,ppn-1):
                          cmd = cmd+'1,'
                      cmd = cmd+'1 --w-cd '+ens_dir
                  elif ( (i % ppn) < (ppn-1)):
                      cmd = cmd+','+ens_dir
                  else:
                    output_run.write(cmd+','+ens_dir+' '+exeroot+'/cesm.exe > ccsm_log.txt &\n')

            if ((i+1) % int(options.ng) == 0 or (i+1) == nsamples):          
                output_copy.write('wait\n')
                output_run.write('wait\n')
                                 
                output_copy.close()
                output_run.close()

                if (not options.no_submit):
                    if (num == 0):
                        if (not options.ensemble_nocopy):
                            subprocess.call('qsub '+tmpdir+'/ensemble_copy_'+ \
                                                casename+'_'+numst[1:]+'.pbs > '+tmpdir+ \
                                                '/jobinfo_copy', shell=True)
                            myinput = open(tmpdir+'/jobinfo_copy')
                            for s in myinput:
                                lastjob = s.split('.')[0]
                            myinput.close()
                            subprocess.call('qsub -W depend=afterok '+lastjob+' '+tmpdir+'/ensemble_run_'+ \
                                                casename+'_'+numst[1:]+'.pbs > '+tmpdir+ \
                                                '/jobinfo_run', shell=True)
                        else:
                            subprocess.call('qsub '+tmpdir+'/ensemble_run_'+ \
                                                casename+'_'+numst[1:]+'.pbs > '+tmpdir+ \
                                                '/jobinfo_run', shell=True) 
                    else:
                        if (not options.ensemble_nocopy):
                            myinput = open(tmpdir+'/jobinfo_copy')
                            for s in myinput:
                                lastjob = s.split('.')[0]
                            myinput.close()
                            subprocess.call('qsub '+tmpdir+'/ensemble_copy_'+ \
                                                casename+'_'+numst[1:]+'.pbs > '+tmpdir+ \
                                                '/jobinfo_copy', shell=True)
                        myinput = open(tmpdir+'/jobinfo_run')
                        for s in myinput:
                            lastjob = s.split('.')[0]
                        myinput.close()
                        subprocess.call('qsub -W depend=afterok:'+lastjob+' '+tmpdir+ \
                              '/ensemble_run_'+casename+'_'+ \
                              numst[1:]+'.pbs > '+tmpdir+'/jobinfo_run', shell=True)
                num = num+1
