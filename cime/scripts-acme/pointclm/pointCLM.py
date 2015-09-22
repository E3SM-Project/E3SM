#!/usr/bin/env python

import os, sys, csv, time, math
from optparse import OptionParser
from Numeric import *


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
parser.add_option("--caseroot", dest="caseroot", default='../../', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--cpl_bypass", dest="cpl_bypass", default=False, \
                   help = "Bypass coupler (point CLM only)", action = "store_true")
parser.add_option("--runroot", dest="runroot", default="../../../run", \
                  help="Directory where the run would be created")
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
parser.add_option("--csmdir", dest="csmdir", default='../../..', \
                  help = "base CESM directory (default = ../../..)")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../../../ccsm_inputdata', \
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
parser.add_option("--npoolmod", action="store_true", dest="npoolmod", default=False, \
                    help="To turn on nitrogen pool modifications")
parser.add_option("--cpoolmod", action="store_true", dest="cpoolmod", default=False, \
                    help="To turn on carbon storage pool modifications")

parser.add_option("--q10wbmod", action="store_true", dest="q10wbmod", default=False, \
                    help="To turn on Woodrow-Berry Q10 curve (CLM 4.0 only)")
parser.add_option("--tfmod", action="store_true", dest="tfmod", default=False, \
                    help="To set temperature threshold (0 degC) for plant wilting factor")
parser.add_option("--harvmod", action="store_true", dest="harvmod", \
                      default=False, help = "Turn on harvest modificaton" \
                      "All harvest is performed in first timestep")
parser.add_option("--humhol", action="store_true", dest="humhol", \
                      default=False, help = "SPRUCE Hummock/Hollow modification")
parser.add_option("--vertsoilc", dest="vsoilc", default=False, \
                  help = 'To turn on CN with multiple soil layers, excluding CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me', action="store_true")
parser.add_option("--MICROBE", dest="MICROBE", default=False, \
                  help = 'To turn on MICROBE with CN', action="store_true")
#parser.add_option("--arcticpft", dest="arcticpft", default=False, \
#                  help = 'To turn on Expanded Arctic PFTs flag (-DPFTARCTIC) in CLM4.5. Must provide --parm_file', action="store_true")
#parser.add_option("--C13", dest="C13", default=False, \
#                  help = 'Switch to turn on C13', action="store_true")
#parser.add_option("--C14", dest="C14", default=False, \
#                  help = 'Use C14 as C13 (no decay)', action="store_true")
parser.add_option("--branch", dest="branch", default=False, \
		  help = 'Switch for branch run', action="store_true")
parser.add_option("--makemetdata", dest="makemet", default=False, \
		  help = 'Generate meteorology', action="store_true")
parser.add_option("--surfdata_grid", dest="surfdata_grid", default=False, \
                  help = 'Use gridded surface data instead of site data', action="store_true")
parser.add_option("--include_nonveg", dest="include_nonveg", default=False, \
                  help = 'Include non-vegetated columns/Landunits in surface data')
parser.add_option("--refcase", dest="refcase" , default='none', \
                  help = 'Use already compiled CLM case')
parser.add_option("--xpts", dest="xpts", default=1, \
                      help = 'for regional runs: xpts')
parser.add_option("--ypts", dest="ypts", default=1, \
                      help = 'for regional runs: ypts')
parser.add_option("--trans2", dest="trans2", default=False, action="store_true", \
                  help = 'Tranisnent phase 2 (1901-2010) - CRUNCEP only')
parser.add_option("--spinup_vars", dest="spinup_vars", default=False, \
                  help = 'Limit output vars in spinup runs', action="store_true")

(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------

PTCLMdir = os.getcwd()

#check for valid csm directory
if (os.path.exists(options.csmdir) == False):
    print('Error:  invalid CESM root directory')
    sys.exit()
else:
    csmdir=os.path.abspath(options.csmdir)
    scriptsdir = csmdir+'/scripts'

#case directory
if (options.caseroot == '' or (os.path.exists(options.caseroot) == False)):
    caseroot = csmdir+'/cases'
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

#check for valid compset
compset = options.compset
#if (compset != 'I1850CLM45CN' and compset != 'I1850CLM45' and compset != 'I2000CLM45' and \
#        compset != 'I2000CLM45CN' and compset != 'I20TRCLM45CN' and compset != 'I1850CN' and \
#        compset != 'I2000CN' and compset != 'I1850' and compset != 'I2000' and \
#        compset != 'I20TRCN' and compset != 'I20TRCRUCLM45CN' and compset !='I1850CLM45BGC' \ 
#        and compset !='I2000CLM45BGC' and compset != 'I20TRCLM45BGC' and \
#        compset != 'I1850CRUCLM45BGC' and compset != 'I20TRCRUCLM45BGC'):
#    print('Error:  Only ICLM40, ICLM45, ICLM40CN, ICLM45CN/BGC compsets are supported')
#    sys.exit()

#figure out if clm40 or clm45, set model-specific options
isclm45 =  False
surfdir = 'surfdata'
pftphys_stamp = 'clm40.c130424'
if ('CLM45' in compset):
    isclm45 = True
    surfdir = 'surfdata_map'
    pftphys_stamp = 'c130821'

#check consistency of options
if ('20TR' in compset):
    #ignore spinup option if transient compset
    if (options.ad_spinup or options.exit_spinup):
      print('Spinup options not available for transient compset.')
      sys.exit()
    #finidat is required for transient compset
    if (options.finidat_case == ''):
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

if (finidat == ''  and options.finidat_case == ''):  #not user-defined
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
    finidat_yst = str(finidat_year)
    if (finidat_year >= 100 and finidat_year < 1000):
        finidat_yst = '0'+str(finidat_year)
    if (finidat_year >= 10 and finidat_year < 100):
        finidat_yst = '00'+str(finidat_year)
    if (finidat_year < 10):
        finidat_yst = '000'+str(finidat_year)
    if(options.runroot == './'): 
        finidat = csmdir+'/run/'+options.finidat_case+'/run/'+ \
                  options.finidat_case+'.clm2.r.'+finidat_yst+ \
                  '-01-01-00000.nc'
    else:
        finidat = runroot+'/'+options.finidat_case+'/run/'+ \
                  options.finidat_case+'.clm2.r.'+finidat_yst+ \
                  '-01-01-00000.nc'

#construct default casename
casename    = options.site+"_"+compset
if (options.mycaseid != ""):
    casename = options.mycaseid+'_'+casename

#CRU-NCEP 2 transient phases
if ('CRU' in compset):
    use_cruncep = True
else:
    use_cruncep = False
if ('20TR' in compset and use_cruncep):
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
blddir=runroot+'/'+casename
print("CASE exeroot is: "+blddir+"\n")
rundir=runroot+'/'+casename+'/run'
print("CASE rundir is: "+rundir+"\n")


#------------------- make point data for site -------------------------------
if (options.nopointdata == False):
    ptcmd = 'python makepointdata.py --caseroot '+options.caseroot+' --casename '+casename+ \
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
os.chdir(csmdir+'/scripts')
#get simyr
mysimyr=1850
if (options.compset == 'ICLM45CN' or options.compset == 'ICLM45BGC' or '2000' in compset):
    mysimyr=2000

#parameter (pft-phys) modifications if desiredi
os.system('mkdir -p '+csmdir+'/scripts/acme/pointclm/temp')
os.system('cp '+options.ccsm_input+'/lnd/clm2/paramdata/clm_params.'+pftphys_stamp+'.nc ' \
              +csmdir+'/scripts/acme/pointclm/temp/clm_params.'+pftphys_stamp+'.'+ \
              casename+'.nc')
os.system('chmod u+w ' +csmdir+'/scripts/acme/pointclm/temp/clm_params.'+pftphys_stamp+'.'+ \
              casename+'.nc')
if (options.parm_file != ''):
    import Scientific.IO.NetCDF
    from Scientific.IO import NetCDF 
    pftfile = NetCDF.NetCDFFile(csmdir+'/scripts/acme/pointclm/temp/' \
                                +'clm_params.'+pftphys_stamp+'.'+casename+'.nc',"a")
    input   = open(os.path.abspath(options.parm_file))
    for s in input:
        if s[0:1] != '#':
            values = s.split()
            temp = pftfile.variables[values[0]]
            temp_data = temp.getValue()
            temp_data[int(values[1])] = float(values[2])
            temp.assignValue(temp_data)
    input.close()
    pftfile.close()

#parameter (soil order dependent) modifications if desired    ::X.YANG 
os.system('cp '+options.ccsm_input+'/lnd/clm2/pftdata/CNP_parameters_c121029.nc ' \
              +csmdir+'/scripts/acme/pointclm/temp/CNP_parameters_c121029'+ \
              casename+'.nc')
os.system('chmod u+w ' +csmdir+'/scripts/acme/pointclm/temp/CNP_parameters_c121029'+ \
              casename+'.nc')
if (options.parm_file_P != ''):
    import Scientific.IO.NetCDF
    from Scientific.IO import NetCDF
    soilorderfile = NetCDF.NetCDFFile(options.ccsm_input+'/lnd/clm2/pftdata/' \
                                +'CNP_parameters_c121029'+casename+'.nc',"a")
    input   = open(os.path.abspath(options.parm_file_P))
    for s in input:
        if s[0:1] != '#':
            values = s.split()
            temp = soilorderfile.variables[values[0]]
            temp_data = temp.getValue()
            temp_data[int(values[1])] = float(values[2])
            temp.assignValue(temp_data)
    input.close()
    soilorderfile.close()

#set number of run years for ad, exit spinup cases
if (options.ny_ad != options.run_n and options.ad_spinup):
    options.run_n = options.ny_ad
elif (options.exit_spinup):
    options.run_n = 1

#------------------IF no refcase, create, configure and build -----------------------------------------

if (options.refcase == 'none'):
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
    #if (options.runroot != ''):
    #    os.system('./xmlchange -file env_build.xml -id EXEROOT -val '+runroot+'/'+casename)

    #turn off ROF module
    os.system('./xmlchange -file env_build.xml -id RTM_MODE -val NULL')

    #clm 4_5 cn config options
    if (isclm45):
        clmcn_opts = "-phys clm4_5"
    else:
        clmcn_opts = "-phys clm4_0"

    #if ("CN" in compset):
    #    clmcn_opts += " -bgc cn"
    #elif ("BGC" in compset):
    #    clmcn_opts += " -bgc bgc"

    #if (options.MICROBE):
    #    clmcn_opts += " -microbe on"

    #if (options.nofire and isclm45):
    #    clmcn_opts += " -nofire"

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
    if (options.runroot != ''):
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
                      +'RUN_REFDATE -val '+finidat_yst+'-01-01')
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_REFCASE -val '+options.finidat_case)
    else:
        if (('CN' in compset or 'BGC' in compset) and options.ad_spinup == False and \
                options.coldstart==False):
            os.system('./xmlchange -file env_run.xml -id RUN_REFDATE -val ' \
                          +finidat_yst+'-01-01')

    #adds capability to run with transient CO2
    if ('20TR' in compset):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'CCSM_BGC -val CO2A')
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'CLM_CO2_TYPE -val diagnostic')

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
            output.write(" hist_mfilt = "+ str(options.hist_mfilt)+"\n")
        if (options.hist_nhtfrq != -999):
            output.write(" hist_nhtfrq = "+ str(options.hist_nhtfrq)+"\n")
        if (options.hist_vars != ''):
            output.write(" hist_empty_htapes = .true.\n")
            #read hist_vars file
            hvars_file = open('../'+options.hist_vars)
            myline = " hist_fincl1 = "
            for s2 in hvars_file:
                if line2 ==0:
                    myline = myline+"'"+s2.strip()+"'"
                else:
                    myline = myline+",'"+s2.strip()+"'"
                line2=line2+1
                output.write(myline+"\n")
                hvars_file.close()
	if (options.spinup_vars):
	  output.write(" hist_empty_htapes = .true.\n")
          output.write(" hist_fincl1 = 'NEE', 'GPP', 'FPSN', 'AR', 'HR', 'MR', 'GR', 'ER', 'NPP', 'TLAI', 'TOTSOMC', 'LEAFC', 'DEADSTEMC', 'DEADCROOTC', 'FROOTC', 'LIVESTEMC', 'LIVECROOTC', 'TOTVEGC', 'TOTCOLC', 'TOTLITC', 'TOTLITN', 'CWDC', 'QVEGE', 'QVEGT', 'QSOIL', 'QDRAI', 'QRUNOFF', 'H2OSFC', 'ZWT', 'FPI', 'FPG'\n")

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
            output.write(" flanduse_timeseries = './surfdata.pftdyn_"+str(numxpts)+'x' \
                 +str(numypts)+"pt_"+casename+".nc'\n")
            output.write(' check_finidat_fsurdat_consistency = .false.\n')
        #pft-physiology file
        output.write(" paramfile = './clm_params."+pftphys_stamp+"."+ \
                     casename+".nc'\n")
        #soil order parameter file
        #output.write(" fsoilordercon = './CNP_parameters_c121029"+ \
        #             casename+".nc'\n")

        #nitrogen deposition file
        if ('CN' in compset or 'BGC' in compset):
            output.write( " stream_fldfilename_ndep = '"+options.ccsm_input+ \
      "/lnd/clm2/ndepdata/fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc'\n")
	if (options.vsoilc):
            output.write(" use_vertsoilc = .true.\n")
	if (options.centbgc):
	    output.write(" use_century_decomp = .true.\n")
	if (options.CH4):
	    output.write(" use_lch4 = .true.\n")
        if (options.nofire and isclm45):
            output.write(" use_nofire = .true.\n")
        output.close()

    #configure case
    if (options.no_config == False):
        os.system('./cesm_setup')
    else:
        print("Warning:  No case configure performed")
        sys.exit()

    #stream file modificaitons: directory and domain file (for using site_level CRU-NCEP)
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
            if (options.npoolmod):
                print("Turning on plant NPOOL modificaiton\n")
	        stemp = stemp[:-1]+' -DNPOOLMOD\n'
            if (options.cpoolmod):
                print("Turning on plant CPOOL modificaiton\n")
	        stemp = stemp[:-1]+' -DCPOOLMOD\n'
            if (options.tfmod):
                print("Turning on TFMOD modification\n")
                stemp = stemp[:-1]+' -DTFMOD\n'
            if (options.q10wbmod):
                print("Turning on Q10WBMOD modification\n")
                stemp = stemp[:-1]+' -DQ10WBMOD\n'
            if (options.harvmod):
                print("Turning on HARVMOD modificaiton\n")
                stemp = stemp[:-1]+' -DHARVMOD\n'
            if (options.humhol):
                print("Turning on HUM_HOL modification\n")
                stemp = stemp[:-1]+' -DHUM_HOL\n'
            if (options.ad_spinup and isclm45 == False):
                print("Turning on AD_SPINUP (CLM4.0)")
                stemp = stemp[:-1]+' -DAD_SPINUP\n'
            if (options.exit_spinup and isclm45 == False): 
                print("Turning on EXIT_SPINUP (CLM4.0)")
                stemp = stemp[:-1]+' -DEXIT_SPINUP\n'
            if (options.cpl_bypass):
                stemp = stemp[:-1]+' -DCPL_BYPASS\n'
                #put the correct input filename in the code
                os.system('cp ../../models/lnd/clm/src/cpl/lnd_import_export.F90 ./SourceMods/src.clm')
                src_input  = open('./SourceMods/src.clm/lnd_import_export.F90','r')
                src_output = open('./SourceMods/src.clm/lnd_import_export_tmp.F90','w')
                for s in src_input:
                    if ('#THISSITEFILE#' in s):
                        ptstr = str(numxpts)+'x'+str(numypts)+'pt'
                        src_output.write("ierr = nf90_open('"+options.ccsm_input+"/atm/datm7/CLM1PT_data/"+ptstr+"_"+options.site+"/all_hourly.nc', NF90_NOWRITE, ncid)\n")
                    elif ('#THISSITEHDMFILE#' in s):
		        src_output.write("ierr = nf90_open('"+options.ccsm_input+"/lnd/clm2/firedata/clmforc.Li_2012_"+ptstr+"_"+options.site+".nc', NF90_NOWRITE, ncid)\n")
	            elif ('#THISSITELNFMFILE#' in s):
		        src_output.write("ierr = nf90_open('"+options.ccsm_input+"/atm/datm7/NASA_LIS/clmforc.Li_2012_climo1995-2011."+ptstr+"_"+options.site+".nc', NF90_NOWRITE, ncid)\n")
                    elif ('#THISSITENDEPFILE#' in s):
                        src_output.write("ierr = nf90_open('"+options.ccsm_input+"/lnd/clm2/ndepdata/fndep_clm_simyr1849-2006_"+ptstr+"_"+options.site+".nc', nf90_nowrite, ncid)\n")
                    elif ('#THISCO2FILE#' in s):
                         src_output.write("ierr = nf90_open('"+options.ccsm_input+"/atm/datm7/CO2/fco2_datm_1765-2007_c100614.nc', nf90_nowrite, ncid)\n")
                    elif ('#THISSITEAEROFILE#' in s):
                         src_output.write("ierr = nf90_open('"+options.ccsm_input+"/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1849-2006_"+ptstr+"_"+options.site+".nc', nf90_nowrite, ncid)\n")
                    elif ('#STARTHOUR#' in s):
                        nyr_cycle = endyear-startyear+1
                        if ('20TR' in compset and options.align_year != -999):
                            starthour = (nyr_cycle - (int(options.align_year) - 1850))*8760
                        else:
                            starthour = nyr_cycle*8760
                        src_output.write(s.replace('#STARTHOUR#',str(starthour)))
                    else:
                        src_output.write(s)
                os.system('mv ./SourceMods/src.clm/lnd_import_export_tmp.F90 ./SourceMods/src.clm/lnd_import_export.F90')
                src_input.close()
                src_output.close()
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
    if(options.caseroot == './' ):
        os.chdir(csmdir+"/scripts/"+casename)
    else:
        os.chdir(casedir)       

    #Datm mods/ transient CO2 patch for transient run (datm buildnml mods)
    myinput  = open('./Buildconf/datmconf/datm_atm_in')
    myoutput = open('user_nl_datm','w')
    for s in myinput:
        if ('streams =' in s):
            myalign_year = 1 #startyear
            if (options.align_year != -999):
                myalign_year = options.align_year
            if ('20TR' in compset):
                mypresaero = '"datm.streams.txt.presaero.trans_1850-2000 1850 1850 2000"'
                myco2      = ', "datm.global1val.streams.co2.txt 1766 1766 2010"'
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
        elif ('taxmode =' in s):
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

    if ('20TR' in compset):  
        os.system('cp '+csmdir+'/models/lnd/clm/doc/UsersGuide/co2_streams.txt ./')
        myinput  = open('co2_streams.txt','r')
        myoutput = open('co2_streams.txt.tmp','w')
        for s in myinput:
            s2 = s.strip()
            if (s2 == '<filePath>'):
                myoutput.write(s)
                myoutput.write('            '+options.ccsm_input+'/atm/datm7/CO2\n')
                next(myinput)
            elif (s2 == '<fileNames>'):
                myoutput.write(s)
                myoutput.write('            '+options.co2_file+'\n')
                next(myinput)
            else:
                myoutput.write(s)
        myinput.close()
        myoutput.close()
        os.system('mv co2_streams.txt.tmp co2_streams.txt')

   #clean build if requested
    if (options.clean_build):
        os.system('./'+casename+'.clean_build')
    #compile cesm
    if (options.no_build == False):
        os.system('./'+casename+'.build')
        if ('20TR' in compset):
            #note:  *.build will sweep everything under Buildconf, but we need 'co2streams.txt'
            #        in transient run
            os.system('cp -f co2_streams.txt ./Buildconf/datmconf/datm.global1val.streams.co2.txt')
            os.system('cp -f co2_streams.txt '+rundir+'/datm.global1val.streams.co2.txt')
    if (options.caseroot == './'):
        os.chdir(csmdir+"/scripts/"+casename)
    else:
        os.chdir(casedir) 


else:  

#----------------------------Reference case set ------------------------------------------
    
    os.chdir('../run')
    incasename  = options.refcase+'_'+options.compset
    if  (use_cruncep and '20TR' in compset):
        if (options.trans2):
            incasename = incasename + '_phase2'
        else:
            incasename = incasename + '_phase1'
    if (options.ad_spinup):
        incasename = incasename + '_ad_spinup'
    if (options.exit_spinup):
        incasename = incasename + '_exit_spinup'
    os.system('mkdir -p '+casename+'/run')
    os.chdir(casename+'/run')
    os.system('pwd')
           
    print 'Copying files from '+incasename+' to '+casename
    os.system('cp '+csmdir+'/run/'+incasename+'/run/*_in* .')
    os.system('cp '+csmdir+'/run/'+incasename+'/bld/cesm.exe .')
    os.system('cp '+csmdir+'/run/'+incasename+'/run/*.nml .')
    os.system('cp '+csmdir+'/run/'+incasename+'/run/*eam* .')
    os.system('cp '+csmdir+'/run/'+incasename+'/run/*.rc .')

   #Change generic site/case name to actual site/case name in namelst files
    os.system('chmod u+w *')
    os.system('sed -e s/'+incasename+'/'+casename+'/ig  lnd_in > lnd_in_tmp')
    os.system('mv lnd_in_tmp lnd_in')
    os.system('sed -e s/REFCASE/'+options.site+'/ig  lnd_in > lnd_in_tmp')
    os.system('mv lnd_in_tmp lnd_in')
    ptstr = str(numxpts)+'x'+str(numypts)+'pt'
    os.system('sed -e s/1x1pt/'+ptstr+'/ig  lnd_in > lnd_in_tmp')
    os.system('mv lnd_in_tmp lnd_in')
    #os.system('sed -e s/'+incasename+'/'+casename+'/ig  datm_atm_in > datm_atm_in_tmp')
    #os.system('mv datm_atm_in_tmp datm_atm_in')
    #os.system('sed -e s/REFCASE/'+options.site+'/ig  datm_atm_in > datm_atm_in_tmp')
    #os.system('mv datm_atm_in_tmp datm_atm_in')
    os.system('sed -e s/1x1pt_REFCASE/'+ptstr+'_'+options.site+'/ig  datm_atm_in > datm_atm_in_tmp')
    os.system('mv datm_atm_in_tmp datm_atm_in')
    #os.system('sed -e s/CLM_USRDAT/1x1pt_'+options.site+'/ig  datm_atm_in > datm_atm_in_tmp')
    #os.system('mv datm_atm_in_tmp datm_atm_in')
    #os.system('mv clm1PT.CLM_USRDAT.stream.txt clm1PT.1x1pt_REFCASE.stream.txt')
    if (use_cruncep):
        os.system('sed -e s/1x1pt_REFCASE/'+ptstr+'_'+options.site+'/ig datm.streams.txt.CLMCRUNCEP.Precip > ' \
                      +'datm.streams.txt.CLMCRUNCEP.Precip.tmp')
        os.system('sed -e s/1x1pt_REFCASE/'+ptstr+'_'+options.site+'/ig datm.streams.txt.CLMCRUNCEP.Solar > ' \
                      +'datm.streams.txt.CLMCRUNCEP.Solar.tmp')
        os.system('sed -e s/1x1pt_REFCASE/'+ptstr+'_'+options.site+'/ig datm.streams.txt.CLMCRUNCEP.TPQW > ' \
                      +'datm.streams.txt.CLMCRUNCEP.TPQW.tmp')
    else:
        os.system('sed -e s/1x1pt_REFCASE/'+ptstr+'_'+options.site+'/ig datm.streams.txt.CLM1PT.CLM_USRDAT > ' \
                      +'datm.streams.txt.CLM1PT.CLM_USRDAT.tmp')
    #os.system('rm *REFCASE*')
    os.system('sed -e s/'+incasename+'/'+casename+'/ig  drv_in > drv_in_tmp')
    os.system('mv drv_in_tmp drv_in')
    os.system('sed -e s/REFCASE/'+options.site+'/ig  drv_in > drv_in_tmp')
    os.system('mv drv_in_tmp drv_in')
    
    #modify met stream file for correct years
    if (use_cruncep):
        types = ['CRUNCEP.Precip', 'CRUNCEP.Solar', 'CRUNCEP.TPQW']
    else:
        types = ['1PT.CLM_USRDAT']
    for t in types:                  
        if (use_cruncep):
            if (t == 'CRUNCEP.Precip'):
                prefix = 'clmforc.cruncep.c2010.0.5d.Prec.'
            elif (t == 'CRUNCEP.Solar'):
                prefix = 'clmforc.cruncep.c2010.0.5d.Solr.'
            elif (t == 'CRUNCEP.TPQW'):
                prefix = 'clmforc.cruncep.c2010.0.5d.TPQWL.'          
        else:
            prefix = ''
        myinput  = open('datm.streams.txt.CLM'+t+'.tmp', 'r')
        myoutput = open('datm.streams.txt.CLM'+t,'w')
        if (use_cruncep and not options.trans2):
            firstfile = '1901-01.nc'
        elif (use_cruncep):
            firstfile = '1921-01.nc'
        else:
            firstfile = '2000-01.nc'
        fieldinfo=False
        for s in myinput:
            if ('fieldInfo' in s):
                fieldinfo = True
                myoutput.write(s)
            elif (firstfile in s):
                print 'TIME INFO', startyear, endyear
                for y in range(startyear,endyear+1):
                    for m in range(1,13):
                        if (m < 10):
                            myoutput.write('            '+prefix+str(y)+'-0'+str(m)+'.nc\n')
                        else:
                            myoutput.write('            '+prefix+str(y)+'-'+str(m)+'.nc\n')
            elif (fieldinfo and '.nc' in s):
                pass
            else:
                myoutput.write(s)
        myinput.close()
        myoutput.close()
        os.system('rm datm.streams.txt.CLM'+t+'.tmp')

    #modify presearo stream file to change to 1850-2000 file
    if ('20TR' not in compset):
        myinput  = open('datm.streams.txt.presaero.clim_1850')
        myoutput = open('datm.streams.txt.presaero.clim_1850.tmp','w')
        for s in myinput:
            if ('aerosoldep_monthly' in s):
                myoutput.write('            aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc\n')
            else:
                myoutput.write(s)
        myinput.close()
        myoutput.close()
        os.system('mv datm.streams.txt.presaero.clim_1850.tmp datm.streams.txt.presaero.clim_1850')

    #modify datm_atm_in for correct years
    myinput  = open('datm_atm_in')
    myoutput = open('datm_atm_in_tmp','w')
    for s in myinput:
        if ('streams =' in s):
            myalign_year = 1 #startyear
            if (options.align_year != -999):
                myalign_year = options.align_year
            if ('20TR' in compset):
                mypresaero = '"datm.streams.txt.presaero.trans_1850-2000 1850 1850 2000"'
                myco2      = ', "datm.global1val.streams.co2.txt 1766 1766 2010"'
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
        elif ('taxmode =' in s):
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
    os.system('mv datm_atm_in_tmp datm_atm_in')

    #modify component .nml files
    nmlfiles=['atm','cpl','glc','ice','lnd','ocn','rof','wav']
    for mynml in nmlfiles:
        outfile = open(mynml+'_modelio.nml','w')
        outfile.write('&modelio\n')
        outfile.write('   diri    = "'+os.path.abspath('../..')+'/'+incasename+'/'+ \
                          mynml+'   "\n')
        outfile.write('   diro    = "./"\n') #'+os.path.abspath('.')+'   "\n')
        outfile.write('   logfile = "'+mynml+'.log   "\n')
        outfile.write('/\n')
        outfile.write('&pio_inparm\n')
        outfile.write(' pio_numiotasks = -99\n')
        outfile.write(' pio_root = -99\n')
        outfile.write(' pio_stride = -99\n')
        outfile.write(" pio_typename = 'nothing'\n")
        outfile.write(" /\n")
        outfile.close()
          
    #make drv_in namelist modifications (run length for final spin/tranisent case)
    myinput  = open('drv_in')
    myoutput = open('drv_in_tmp','w')
    for s in myinput:
        if (s[0:7] == ' stop_n'):
            myoutput.write(" stop_n = "+str(options.run_n)+'\n')                               
        elif (s[0:10] == ' restart_n'):
            myoutput.write(" restart_n = "+str(options.run_n)+'\n')
        elif (s[0:9] == ' stop_ymd'):
            myoutput.write(" stop_ymd       = -999\n")
        elif (s[0:12] == ' restart_ymd'):
            myoutput.write(" restart_ymd    = -999\n")
        elif (s[0:11] == ' atm_cpl_dt'):
            myoutput.write(" atm_cpl_dt = "+str(int(float(options.tstep)*3600.0))+'\n')
        elif (s[0:11] == ' lnd_cpl_dt'):
            myoutput.write(" lnd_cpl_dt = "+str(int(float(options.tstep)*3600.0))+'\n')
        elif (s[0:11] == ' ice_cpl_dt'):
            myoutput.write(" atm_cpl_dt = "+str(int(float(options.tstep)*3600.0))+'\n')
        elif (s[0:10] == ' start_ymd'):
            myoutput.write(s)
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()
    os.system('mv drv_in_tmp drv_in')
            
    #make lnd_in namelist modification for finidat file in transient case for correct year
    myinput  = open('lnd_in')
    myoutput = open('lnd_in_tmp','w')
   
    for s in myinput:
        if (s[0:8] == ' hist_mf'):
            if ('20TR' in compset):
                myoutput.write(' hist_mfilt = 12, 8760, 365\n')
            else: 
                myoutput.write(' hist_mfilt = 1\n')
        elif (s[0:8] == ' hist_nh'):
            if ('20TR' in compset):
                myoutput.write(' hist_nhtfrq = 0, -1, -24\n')
                myoutput.write(" hist_fincl2 = 'NEE', 'GPP', 'NPP', 'ER', 'AR', 'MR', 'GR', 'SR', 'EFLX_LH_TOT', 'FSH', 'FPSN', 'BTRAN', 'FPG', 'FPI', 'CPOOL', 'NPOOL', 'TV', 'FSA', 'FIRA', 'FCTR', 'FCEV', 'FGEV', 'TBOT', 'FLDS', 'FSDS', 'RAIN', 'SNOW', 'WIND', 'PBOT', 'QBOT'\n")
                myoutput.write(" hist_fincl3 = 'NEE', 'GPP', 'NPP', 'AGNPP', 'BGNPP', 'ER', 'AR', 'MR', 'GR', 'HR', 'SR', 'EFLX_LH_TOT', 'FSH', 'FPSN', 'BTRAN', 'FPG', 'FPI', 'TBOT', 'FLDS', 'FSDS', 'RAIN', 'SNOW', 'WIND', 'PBOT', 'QBOT', 'PFT_FIRE_CLOSS', 'LITFALL', 'TLAI', 'LEAFC' ,'FROOTC', 'LIVESTEMC', 'DEADSTEMC', 'LIVECROOTC', 'DEADCROOTC', 'TOTVEGC', 'TOTSOMC', 'TOTLITC', 'CWDC', 'TOTECOSYSC', 'TOTCOLC', 'TOTSOMN', 'TOTECOSYSN', 'SMINN', 'QOVER', 'QDRAI', 'QRGWL', 'QRUNOFF' \n")
            else:
                myoutput.write(' hist_nhtfrq = -8760\n')
        elif (s[0:8] == ' finidat' and options.ad_spinup == False):
            myoutput.write(" finidat = '"+finidat+"'\n")
        elif (s[0:7] == ' nrevsn'):
            myoutput.write(" nrevsn = '"+finidat+"'\n")
        elif (s[0:6] == ' dtime'):
            myoutput.write(" dtime = "+str(int(float(options.tstep)*3600.0))+'\n')
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()
    os.system('mv lnd_in_tmp lnd_in')

    #write a basic PBS script
    output = open(casename+'.run','w')
    output.write('#PBS -S /bin/bash\n')
    output.write('#PBS -V\n')
    output.write('#PBS -m ae\n')
    output.write('#PBS -N '+casename+'\n')
    output.write('#PBS -q esd08q\n')
    output.write("#PBS -l nodes="+str((int(options.np)-1)/8+1)+ \
                     ":ppn="+str(min(int(options.np),8))+"\n")  
    output.write('#PBS -l walltime=48:00:00\n')
    output.write("cd "+csmdir+'/run/'+casename+"/run\n")
    if (options.np == 1):
        output.write("./cesm.exe > cesm_log.txt\n")
    else:
        output.write("mpirun -np "+options.np+" --hostfile $PBS_NODEFILE ./cesm.exe\n")
    output.close()

#---------------------------end of refcase ------------------------------------------------



#copy rpointers and restart files to current run directory
if(finidat != '' and options.runroot != '' ):
    os.system('cp -f '+finidat+' '+runroot+'/'+\
              casename+'/run/')
    os.system('cp -f '+runroot+'/'+options.finidat_case+\
              '/run/rpointer.* '+options.runroot+'/run/'+
              casename+'/run')        
if (finidat != '' and options.runroot == '' ):
    os.system('cp -f '+csmdir+'/run/'+options.finidat_case+'/run/'+ \
              options.finidat_case+'.*'+finidat_yst+'* '+csmdir+ \
              '/run/' +casename+'/run/')
    os.system('cp -f '+csmdir+'/run/'+options.finidat_case+'/run/'+ \
              'rpointer.* '+csmdir+'/run/'+casename+'/run/')
#move site data to run directory
os.system('mv '+csmdir+'/scripts/acme/pointclm/temp/*'+casename+'* ' \
           +runroot+'/'+casename+'/run/')
if (options.nopointdata == False):
   os.system('mv '+csmdir+'/scripts/acme/pointclm/temp/domain*'+options.site+'* ' \
          +runroot+'/'+casename+'/run/')


#os.system('cp -f ../microbepar_in ' +csmdir+'/run/'+casename+'/run/')

#submit job if requested
if (options.no_submit == False):
    os.system("pwd")
    os.system("qsub "+casename+".run")


#copy call_PTCLM.py to case directory
#os.chdir('..')
#os.system("cp "+cmd+" ./"+casename+'/call_PTCLM_'+casename+'.cmd')
