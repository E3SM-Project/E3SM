#!/usr/bin/env python3

import getpass, os, sys, csv, math
from optparse import OptionParser
import subprocess
import numpy
import re

parser = OptionParser();

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='', \
                  help = "input data directory for CESM (required)")
parser.add_option("--exeroot", dest="exeroot", default="", \
                 help="Location of executable")
parser.add_option("--model_root", dest="csmdir", default='', \
                  help = "base CESM directory")
parser.add_option("--clean_build", action="store_true", default=False, \
                  help="Perform a clean build")
parser.add_option("--compiler", dest="compiler", default = '', \
                  help = "compiler to use (pgi*, gnu)")
parser.add_option("--debugq", dest="debug", default=False, \
                 action="store_true", help='Use debug queue and options')
parser.add_option("--dailyrunoff", dest="dailyrunoff", default=False, \
                 action="store_true", help="Write daily output for hydrology")
parser.add_option("--ilambvars", dest="ilambvars", default=False, \
                 action="store_true", help="Write special outputs for diagnostics")
parser.add_option("--dailyvars", dest="dailyvars", default=False, \
                 action="store_true", help="Write daily ouptut variables")
parser.add_option("--hist_mfilt_trans", dest="hist_mfilt", default="1", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_mfilt_spinup", dest="hist_mfilt_spinup", default="-999", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_nhtfrq_spinup", dest="hist_nhtfrq_spinup", default="-999", \
                  help = 'output file timestep (transient only)')
parser.add_option("--hist_nhtfrq_trans", dest="hist_nhtfrq", default="0", \
                  help = 'output file timestep (transient only)')
parser.add_option("--point_list", dest="point_list", default='', \
                  help = 'File containing list of points to run')
parser.add_option("--lat_bounds", dest="lat_bounds", default='-90,90', \
                  help = 'latitude range for regional run')
parser.add_option("--lon_bounds", dest="lon_bounds", default='-180,180', \
                  help = 'longitude range for regional run')
parser.add_option("--mask", dest="mymask", default='', \
                  help = 'Mask file to use (regional only)')
parser.add_option("--machine", dest="machine", default = '', \
                  help = "machine to use")
parser.add_option("--monthly_metdata", dest="monthly_metdata", default = '', \
                  help = "File containing met data")
parser.add_option("--mpilib", dest="mpilib", default="", \
                      help = "mpi library (openmpi*, mpich, ibm, mpi-serial)")
parser.add_option("--noad", action="store_true", dest="noad", default=False, \
                  help='Do not perform ad spinup simulation')
parser.add_option("--nofire", dest="nofire", default=False, \
                  action="store_true", help='Turn off fire algorithms')
parser.add_option("--nofn", action="store_true", dest="nofn", default=False, \
                  help="Do not perform final spinup simulation")
parser.add_option("--notrans", action="store_true", dest="notrans", default=False, \
                  help='Do not perform transient simulation (spinup only)')
parser.add_option("--nopointdata", action="store_true", dest="nopointdata", \
                  default=False, help="do not generate point data")
parser.add_option("--nyears_final_spinup", dest="nyears_final_spinup", default='200', \
                  help="base no. of years for final spinup")
parser.add_option('--project', dest='project',default='', \
                 help='Set project')
parser.add_option('--region', dest="region", default='', \
                 help="Set pre-defined region")
parser.add_option('--runblock', dest="runblock", default=9999, \
                  help="Number of years to run for each submission")
parser.add_option("--runroot", dest="runroot", default="", \
                  help="Directory where the run would be created")
parser.add_option("--run_startyear", dest="run_startyear",default=-1, \
                      help='Starting year for model output (SP only)')
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--surffile", dest="surffile", default='', \
                  help = 'Use specified surface data file')
parser.add_option("--parm_file", dest="parm_file", default="", \
                  help = 'parameter file to use')
parser.add_option("--parm_file_P", dest="parm_file_P", default="", \
                  help = 'P parameter file to use')               
parser.add_option("--parm_vals", dest="parm_vals", default="", \
                  help = 'User specified parameter values')
parser.add_option("--pft", dest="mypft", default=-1, \
                  help = 'Use this PFT for all gridcells')
parser.add_option("--nopftdyn", action="store_true", dest="nopftdyn", \
                      default=False, help='Do not use dynamic PFT file')
parser.add_option("--res", dest="res", default="hcru_hcru", \
                      help='Resoultion for global simulation')
parser.add_option("--np", dest="np", default=256, \
                  help = 'number of processors')
parser.add_option("--tstep", dest="tstep", default=0.5, \
                  help = 'CLM timestep (hours)')
parser.add_option("--co2_file", dest="co2_file", default="fco2_datm_1765-2007_c100614.nc", \
                  help = 'CLM timestep (hours)')
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=250, \
                  help = 'number of years to run ad_spinup')
parser.add_option("--nyears_transient", dest="nyears_transient", default=-1, \
                  help = 'number of years to run transient')
parser.add_option("--metdir", dest="metdir", default="none", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--C13", dest="C13", default=False, \
                      help = 'Switch to turn on C13', action="store_true")
parser.add_option("--C14", dest="C14", default=False, \
                  help = 'Use C14 as C13 (no decay)', action="store_true")
parser.add_option("--BGC", action="store_true", dest="bgc", default=False, \
                  help = 'Use BGC compset')
parser.add_option("--SP", dest="sp", default=False, \
                 action="store_true", help="Use satellite phenology")
parser.add_option("--harvmod", action="store_true", dest='harvmod', default=False, \
                    help="turn on harvest modification:  All harvest at first timestep")
parser.add_option("--bulk_denitrif", dest="bulk_denitrif", default=False, \
                  help = 'To turn on BGC nitrification-denitrification', action="store_true")
parser.add_option("--no_dynroot", dest="no_dynroot", default=False, \
                  help = 'Turn off dynamic root distribution', action="store_true")
parser.add_option("--vertsoilc", dest="vsoilc", default=False, \
                  help = 'To turn on CN with multiple soil layers, excluding CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--ECA", action="store_true", dest="eca", default=False, \
                  help = 'Use ECA compset')
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me', action="store_true")
parser.add_option("--cruncep", dest="cruncep", default=False, \
                  action="store_true", help = 'Use CRU-NCEP meteorology')
parser.add_option("--cruncepv8", dest="cruncepv8", default=False, \
                  help = "use cru-ncep data", action="store_true")
parser.add_option("--cplhist", dest="cplhist", default=False, \
                  help= "use CPLHIST forcing", action="store_true")
parser.add_option("--site_forcing", dest="site_forcing", default="", \
                   help="Use forcing from site for all gridcells")
parser.add_option("--gswp3", dest="gswp3", default=False, \
                  action="store_true", help = 'Use GSWP3 meteorology')
parser.add_option("--princeton", dest="princeton", default=False, \
                  action="store_true", help = 'Use Princeton meteorology')
parser.add_option("--livneh", dest="livneh", default=False, \
                  action="store_true", help = "Livneh correction to CRU precip (CONUS only)")
parser.add_option("--daymet", dest="daymet", default=False, \
                  action="store_true", help = "Daymet correction to GSWP3 precip (CONUS only)")
parser.add_option("--cpl_bypass", dest = "cpl_bypass", default=False, \
                  help = "Bypass coupler", action="store_true")
parser.add_option("--spinup_vars", dest = "spinup_vars", default=False, \
                  help = "limit output variables for spinup", action="store_true")
parser.add_option("--trans_varlist", dest = "trans_varlist", default='', help = "Transient outputs")
parser.add_option("--cn_only", dest="cn_only", default=False, \
                  help='Carbon/Nitrogen only (saturated P)', action ="store_true")
#Changed by Ming for mesabi
parser.add_option("--archiveroot", dest="archiveroot", default='', \
                  help = "archive root directory only for mesabi")
#Added by Kirk to include the modified parameter file
parser.add_option("--mod_parm_file", dest="mod_parm_file", default='', \
                  help = "adding the path to the modified parameter file")
parser.add_option("--mod_parm_file_P", dest="mod_parm_file_P", default='', \
                  help = "adding the path to the modified parameter file")
parser.add_option("--walltime", dest="walltime", default=24, \
                  help = "desired walltime for each job (hours)")

(options, args) = parser.parse_args()

#------------ define function for pbs submission
def submit(fname, submit_type='qsub', job_depend=''):
    if (job_depend != ''):
        if (submit_type == 'qsub'):
            os.system(submit_type+' -W depend=afterok:'+job_depend+' '+fname+' > temp/jobinfo')
        elif (submit_type == 'sbatch'):
            os.system(submit_type+' -d afterok:'+job_depend+' '+fname+' > temp/jobinfo')
    else:
        os.system(submit_type+' '+fname+' > temp/jobinfo')
    myinput = open('temp/jobinfo')
    for s in myinput:
        thisjob = re.search('[0-9]+', s).group(0)
    myinput.close()
    os.system('rm temp/jobinfo')
    return thisjob

def get_regional_bounds(myregion):
    if (myregion == 'noam'):  #North America
        bounds = [-170.25,-45.25,9.75,79.75]
    elif (myregion == 'bona'):  #Boreal North America
        bounds = [-170.25,-60.25,49.75,79.75]
    elif (myregion == 'tena'):  #Temperate North America
        bounds = [-125.25,-66.25,30.25,49.75]
    elif (myregion == 'conus'):  #CONUS
        bounds = [-125.25,-66.25,23.25,54.75]
    elif (myregion == 'columbia'):   #Columbia river watershed
        bounds = [-126, -108, 40.0, 55.0]
    elif (myregion == 'ceam'):  #Central America
        bounds = [-115.25,-80.25,9.75,30.25]
    elif (myregion == 'soam'):    #South America
        bounds = [-80.25,-40.25,-59.75,12.75]
    elif (myregion == 'nhsa'):  #Northern hemisphere South America
        bounds = [-80.25,-50.25,0.25,12.75]
    elif (myregion == 'shsa'):  #Southern hemisphere South America
        bounds = [-80.25,-40.25,-59.75,0.25]
    elif (myregion == 'euro'):  #Europe
        bounds = [-10.25,30.25,35.25,70.25]
    elif (myregion == 'mide'):  #Middle East
        bounds = [-10.25,60.25,20.24,40.25]
    elif (myregion == 'afrc'):    #Africa
        bounds = [-20.25,45.25,-34.75,20.25]
    elif (myregion == 'nhaf'):  #Northern hemisphere Africa
        bounds = [-20.25,45.25,0.25,20.25]
    elif (myregion == 'shaf'):  #Southern hemisphere Africa
        bounds = [10.25,45.25,-34.75,0.25]
    elif (myregion == 'asia'):  #Asia
        bounds = [30.25,179.75,-10.25,70.25]
    elif (myregion == 'boas'):  #Boreal asia
        bounds = [30.25,179.25,54.75,70.25]
    elif (myregion == 'ceas'):  #Central Asia
        bounds = [30.25,142.58,30.25,54.75]
    elif (myregion == 'seas'):  #Southeast Asia
        bounds = [65.25,120.25,5.25,30.25]
    elif (myregion == 'eqas'):  #Equatorial Asia
        bounds = [99.75,150.25,-10.25,10.25]
    elif (myregion == 'aust'):  #Australia
        bounds = [112.00,154.00,-41.25,-10.50]
    else:
        print('Region not recognized.  Setting bounds to global')
        bounds = [-180,180,-90,90]   #global
    return bounds

#ILAMB diagnostic variables and other required variables for analysis
ilamb_outputs = ['FAREA_BURNED', 'CWDC', 'LEAFC', 'TOTLITC', 'STORVEGC', 'LIVESTEMC', 'DEADSTEMC', \
                 'TOTPRODC', 'FROOTC', 'LIVECROOTC', 'DEADCROOTC', 'SOIL1C', 'SOIL2C', 'SOIL3C', \
                 'SOIL4C', 'TOTSOMC', 'TOTVEGC', 'WOODC', 'QSOIL', 'QVEGE', 'COL_FIRE_CLOSS', \
                 'LITR1C_TO_SOIL1C', 'LITR2C_TO_SOIL2C', 'LITR3C_TO_SOIL3C', 'LAND_USE_FLUX', \
                 'LITFALL', 'GPP', 'FGR', 'TLAI', 'SNOWLIQ', 'SOILICE', 'SOILLIQ', 'QRUNOFF', \
                 'QOVER', 'SOILWATER_10CM', 'NBP', 'LEAFC_ALLOC', 'WOODC_ALLOC', 'QINTR', \
                 'AR', 'GR', 'HR', 'MR', 'FSNO', 'SNOWDP', 'QMELT', 'H2OSNO', 'SNOWBCMSL', \
                 'SNODSTMSL', 'SNOOCMSL', 'QVEGT', 'TSOI', 'WIND', 'EFLX_LH_TOT', 'FCTR', \
                 'FCEV', 'FGEV', 'FSH', 'RH2M', 'Q2M', 'RAIN', 'SNOW', 'PBOT', 'FLDS', 'FIRE', \
                 'FSDS', 'FSR', 'TSA', 'CPOOL', 'NPOOL','FPI','FPG','FPI_P','FPG_P','TOTSOMC_1m']

#----------------------------------------------------------

if (options.ccsm_input != ''):
    ccsm_input = options.ccsm_input
elif (options.machine == 'titan' or options.machine == 'eos'):
    ccsm_input = '/lustre/atlas/world-shared/cli900/cesm/inputdata'
elif (options.machine == 'cades' or options.machine == 'metis'):
    ccsm_input = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/ACME_inputdata/'
elif (options.machine == 'edison' or 'cori' in options.machine):
    ccsm_input = '/project/projectdirs/acme/inputdata'
elif ('anvil' in options.machine):
    ccsm_input = '/home/ccsm-data/inputdata'
elif ('compy' in options.machine):
    ccsm_input = '/compyfs/inputdata'

print options.machine
#default compilers
if (options.compiler == ''):
    if (options.machine == 'titan' or options.machine == 'metis'):
        options.compiler = 'pgi'
    if (options.machine == 'eos' or options.machine == 'edison' or 'cori' in options.machine):
        options.compiler = 'intel'
    if (options.machine == 'cades'):
        options.compiler = 'gnu'
    if (options.machine == 'compy'): 
        options.compiler = 'intel'

#default MPIlibs
if (options.mpilib == ''):    
    if ('cori' in options.machine or 'edison' in options.machine):
        options.mpilib = 'mpt'  
    elif ('cades' in options.machine):
        options.mpilib = 'openmpi'
    elif ('anvil' in options.machine):
        options.mpilib = 'mvapich'
    elif ('compy' in options.machine):
        options.mpilib = 'impi'

print options.mpilib
mycaseid   = options.mycaseid
srcmods    = options.srcmods_loc

translen = int(options.nyears_transient)

if (options.bgc):
    mybgc='BGC'
else:
    mybgc='CN'

#csmdir = os.getcwd()+'/'+options.csmdir
csmdir = options.csmdir

#get project information
myuser = getpass.getuser()
if (options.project != ''):
  myproject = options.project
else: 
  if (os.path.exists(os.environ.get('HOME')+'/.cesm_proj')):
    print 'Getting project from '+os.environ.get('HOME')+'/.cesm_proj'
    myinput = open(os.environ.get('HOME')+'/.cesm_proj','r')
    for s in myinput:
        myproject=s[:-1]
    print 'Project = '+myproject
 
#case run and case root directories
if (options.runroot == ''):
    if (options.machine == 'titan' or options.machine == 'eos'):
        runroot='/lustre/atlas/scratch/'+myuser+'/'+myproject
    elif (options.machine == 'cades' or options.machine == 'metis'):
        runroot='/lustre/or-hydra/cades-ccsi/scratch/'+myuser
    elif (options.project == '' and 'cori' in options.machine or 'edison' in options.machine):
        runroot=os.environ.get('CSCRATCH')+'/acme_scratch/'+options.machine+'/'
    elif ('anvil' in options.machine):
        runroot="/lcrc/group/acme/"+myuser
    elif ('compy' in options.machine):
        runroot='/compyfs/'+myuser+'/e3sm_scratch'
    else:
        runroot = csmdir+'/run'
else:
    runroot = os.path.abspath(options.runroot)
if (options.caseroot == '' or (os.path.exists(options.caseroot) == False)):
    caseroot = os.path.abspath(csmdir+'/cime/scripts')
else:
    caseroot = os.path.abspath(options.caseroot)

if (options.cruncep or options.cruncepv8 or options.gswp3 or options.princeton):
    startyear = 1901
    endyear = 1920
    if (options.cruncep):
       site_endyear = 2013
    elif (options.cruncepv8):
       site_endyear = 2016   
    elif (options.gswp3):
       site_endyear = 2010
    elif (options.princeton):
       site_endyear = 2012
    if (options.livneh):
        startyear = 1950
        endyear = 1969
    if (options.daymet):
        startyear = 1980
        endyear = 2010
elif (options.site_forcing):
   #UMB only - test case
   startyear=2000
   endyear=2014
   site_endyear=2014
else:
    #Qian input data
    startyear = 1948
    endyear   = 1972
    site_endyear = 2004

ncycle   = endyear-startyear+1   #number of years in met cycle
ny_ad = options.ny_ad
ny_fin = options.nyears_final_spinup
if (int(options.ny_ad) % ncycle != 0):
    #AD spinup and final spinup lengths must be multiples of met data cyle.
    ny_ad = str(int(ny_ad) + ncycle - (int(ny_ad) % ncycle))
if (int(options.nyears_final_spinup) % ncycle !=0):
    ny_fin = str(int(ny_fin) + ncycle - (int(ny_fin) % ncycle))

if (translen == -1):
    translen = endyear-1850+1        #length of transient run
    if (options.cpl_bypass and (options.cruncep or options.cruncepv8 or options.gswp3 or options.princeton)):
        translen = site_endyear-1850+1

fsplen = int(ny_fin)
if (options.sp):
    #no spinup, just run over 
    fsplen = site_endyear-startyear+1
 
#get align_year
year_align = (endyear-1850+1) % ncycle

#get regional information
isregional = False
if (options.region != ''):
    mybounds = get_regional_bounds(options.region)
    lon_bounds = mybounds[0:2]
    lat_bounds = mybounds[2:4]
else:
    lat_bounds = options.lat_bounds.split(',')
    lon_bounds = options.lon_bounds.split(',')
lat_bounds = [float(l) for l in lat_bounds]
lon_bounds = [float(l) for l in lon_bounds]
if (lon_bounds[0] > -180 or lon_bounds[1] < 180 or lat_bounds[0] > -90 or \
        lat_bounds[1] < 90):
    isregional=True

basecmd = 'python runcase.py --surfdata_grid --ccsm_input '+ \
    os.path.abspath(ccsm_input)+' --rmold --no_submit --machine ' \
    +options.machine+' --res '+options.res+' --model_root '+csmdir
if (options.point_list != ''):
    basecmd = basecmd+' --point_list '+options.point_list
    lat_bounds = [-999,-999]
    lon_bounds = [-999,-999]
if (srcmods != ''):
    srcmods    = os.path.abspath(srcmods)
    basecmd = basecmd+' --srcmods_loc '+srcmods
if (mycaseid != ''):
    basecmd = basecmd+' --caseidprefix '+mycaseid
if (options.nopointdata):
    basecmd = basecmd+' --nopointdata'
if (options.project != ''):
    basecmd = basecmd+' --project '+options.project
if (options.parm_vals != ''):
    basecmd = basecmd+' --parm_vals '+options.parm_vals
if (options.clean_build):
    basecmd = basecmd+' --clean_build '
if (options.metdir !='none'):
    basecmd = basecmd+' --metdir '+options.metdir
#if (not isregional):
#    print 'TEST'
#    basecmd = basecmd+' --nopointdata '
if (options.C13):
    basecmd = basecmd+' --C13 '
if (options.C14):
    basecmd = basecmd+' --C14 '
if (options.debug):
    basecmd = basecmd+' --debugq'
if (options.monthly_metdata != ''):
    basecmd = basecmd+' --monthly_metdata '+options.monthly_metdata
if (options.nofire):
    basecmd = basecmd+' --nofire'
if (options.harvmod):
    basecmd = basecmd+' --harvmod'
if (int(options.mypft) >= 0):
    basecmd = basecmd+' --pft '+str(options.mypft)
if (options.nopftdyn):
    basecmd = basecmd+' --nopftdyn'
if (options.no_dynroot):
    basecmd = basecmd+' --no_dynroot'
if (options.bulk_denitrif):
    basecmd = basecmd+' --bulk_denitrif'
if (options.vsoilc):
    basecmd = basecmd+' --vertsoilc'
if (options.centbgc):
    basecmd = basecmd+' --centbgc'
if (options.cn_only):
    basecmd = basecmd+' --cn_only'
if (options.CH4):
    basecmd = basecmd+' --CH4'
if (options.site_forcing != ''):
    basecmd = basecmd+' --site_forcing '+options.site_forcing
if (options.cruncep):
    basecmd = basecmd+' --cruncep'
if (options.cruncepv8):
    basecmd = basecmd+' --cruncepv8'
if (options.gswp3):
    basecmd = basecmd+' --gswp3'
if (options.princeton):
    basecmd = basecmd+' --princeton'
if (options.livneh):
    basecmd = basecmd+' --livneh'
if (options.daymet):
    basecmd = basecmd+' --daymet'
if (options.cplhist):
    basecmd = basecmd+' --cplhist'
if (options.mymask != ''):
    basecmd = basecmd+' --mask '+options.mymask
if (options.archiveroot !=''):
    basecmd = basecmd+' --archiveroot '+options.archiveroot
if (options.parm_file !=''):
    basecmd = basecmd+' --parm_file '+options.parm_file
if (options.parm_file_P !=''):
    basecmd = basecmd+' --parm_file_P '+options.parm_file_P
if (options.mod_parm_file !=''):
    basecmd = basecmd+' --mod_parm_file '+options.mod_parm_file
if (options.mod_parm_file_P !=''):
    basecmd = basecmd+' --mod_parm_file_P '+options.mod_parm_file_P
if (options.surffile != ''):
    basecmd = basecmd+' --surffile '+options.surffile
basecmd = basecmd + ' --np '+str(options.np)
basecmd = basecmd + ' --tstep '+str(options.tstep)
basecmd = basecmd + ' --co2_file '+options.co2_file
if (options.compiler != ''):
    basecmd = basecmd + ' --compiler '+options.compiler
if (options.mpilib != ''):
    basecmd = basecmd + ' --mpilib '+options.mpilib
basecmd = basecmd+' --caseroot '+caseroot
basecmd = basecmd+' --runroot '+runroot
basecmd = basecmd+' --walltime '+str(options.walltime)
basecmd = basecmd+' --lat_bounds '+str(lat_bounds[0])+','+str(lat_bounds[1])
basecmd = basecmd+' --lon_bounds '+str(lon_bounds[0])+','+str(lon_bounds[1])

#----------------------- build commands for runCLM.py -----------------------------

#ECA or CTC
if (options.cn_only):
    nutrients = 'CN'
else:
    nutrients = 'CNP'
if (options.centbgc):
    decomp_model = 'CNT'
else:
    decomp_model = 'CTC'
if (options.eca):
    mymodel = nutrients+'ECA'+decomp_model
else:
    mymodel = nutrients+'RD'+decomp_model
if (options.cpl_bypass):
    compset_type = 'ICB'
else:
    compset_type = 'I'
mymodel_fnsp = compset_type+'1850'+mymodel+'BC'
mymodel_adsp = mymodel_fnsp.replace('CNP','CN')
mymodel_trns = mymodel_fnsp.replace('1850','20TR')
if (options.sp):
    mymodel_fnsp = compset_type+'CLM45BC'
    options.noad = True
    options.notrans = True

#AD spinup
res=options.res

cmd_adsp = basecmd+' --ad_spinup --nyears_ad_spinup '+ \
    str(ny_ad)+' --align_year '+str(year_align+1)
if (int(options.hist_mfilt_spinup) == -999):
    cmd_adsp = cmd_adsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
        str((endyear-startyear+1)*8760)+' --compset '+mymodel_adsp
else:
    cmd_adsp = cmd_adsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
        +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)+' --compset '+ \
        mymodel_adsp
if (options.exeroot != ''):
  ad_exeroot = os.path.abspath(options.exeroot)
  cmd_adsp = cmd_adsp+' --no_build --exeroot '+ad_exeroot
else:
  ad_exeroot = os.path.abspath(runroot+'/'+ad_case+'/bld')

ad_case = res+'_'+mymodel_adsp+'_ad_spinup'

if (options.spinup_vars):
    cmd_adsp = cmd_adsp+' --spinup_vars'
if (mycaseid != ''):
    ad_case = mycaseid+'_'+ad_case


#final spinup
if mycaseid !='':
    basecase=mycaseid+'_'+res
    basecase = basecase+'_'+mymodel_fnsp
else:
    basecase=res+'_'+mymodel_fnsp
if (options.noad):
    if (options.sp):
        if (options.run_startyear < 0):
            options.run_startyear = startyear
        cmd_fnsp = basecmd+' --run_units nyears --run_n '+str(fsplen)+' --align_year '+ \
            str(year_align+1)+' --coldstart --run_startyear '+str(options.run_startyear)
    else:
        cmd_fnsp = basecmd+' --run_units nyears --run_n '+str(fsplen)+' --align_year '+ \
            str(year_align+1)+' --coldstart'
else:
    cmd_fnsp = basecmd+' --finidat_case '+ad_case+' '+ \
        '--finidat_year '+str(int(ny_ad)+1)+' --run_units nyears --run_n '+ \
        str(fsplen)+' --align_year '+str(year_align+1)+' --no_build' + \
        ' --exeroot '+ad_exeroot+' --nopointdata'
if (int(options.hist_mfilt_spinup) == -999):
    cmd_fnsp = cmd_fnsp+' --hist_mfilt 1 --hist_nhtfrq -'+ \
        str((endyear-startyear+1)*8760)+' --compset '+mymodel_fnsp
else:
    cmd_fnsp = cmd_fnsp+' --hist_mfilt '+str(options.hist_mfilt_spinup) \
        +' --hist_nhtfrq '+str(options.hist_nhtfrq_spinup)+' --compset ' \
        +mymodel_fnsp
if (options.spinup_vars):
    cmd_fnsp = cmd_fnsp+' --spinup_vars'
if (options.dailyrunoff):
    cmd_fnsp = cmd_fnsp+' --dailyrunoff'

#transient
cmd_trns = basecmd+' --finidat_case '+basecase+ \
    ' --finidat_year '+str(fsplen+1)+' --run_units nyears' \
    +' --run_n '+str(translen)+' --align_year '+ \
    str(year_align+1850)+' --hist_nhtfrq '+ \
    options.hist_nhtfrq+' --hist_mfilt '+options.hist_mfilt+' --no_build' + \
    ' --exeroot '+ad_exeroot+' --nopointdata --compset '+mymodel_trns
if (options.trans_varlist != ''):
    cmd_trns = cmd_trns + ' --trans_varlist '+options.trans_varlist
if (options.ilambvars):
    cmd_trns = cmd_trns + ' --ilambvars'
if (options.dailyvars):
    cmd_trns = cmd_trns + ' --dailyvars'        
if (options.dailyrunoff):
    cmd_trns = cmd_trns+' --dailyrunoff'

#transient phase 2 (CRU-NCEP only, without coupler bypass)
if ((options.cruncep or options.gswp3 or options.cruncepv8) and not options.cpl_bypass):
    basecase=basecase.replace('1850','20TR')+'_phase1'
    thistranslen = site_endyear - 1921 + 1
    cmd_trns2 = basecmd+' --trans2 --finidat_case '+basecase+ \
        ' --finidat_year 1921 --run_units nyears --branch ' \
        +' --run_n '+str(thistranslen)+' --align_year 1921'+ \
        ' --hist_nhtfrq '+options.hist_nhtfrq+' --hist_mfilt '+ \
        options.hist_mfilt+' --no_build'+' --exeroot '+ad_exeroot + \
        ' --compset '+mymodel_trns+' --nopointdata'
    print(cmd_trns2)

#---------------------------------------------------------------------------------

os.system('mkdir -p temp')

#build cases
if (options.noad == False):
    print('\nSetting up ad_spinup case\n')
    os.system(cmd_adsp)
if (options.nofn == False):
    print('\nSetting up final spinup case\n')
    os.system(cmd_fnsp)
if (options.notrans == False):
    print('\nSetting up transient case\n')
    os.system(cmd_trns)
if ((options.cruncep or options.gswp3 or options.cruncepv8) and not options.cpl_bypass):
    print('\nSetting up transient case phase 2\n')
    os.system(cmd_trns2)
        

cases=[]
job_depend_run=''
if mycaseid !='':
    basecase=mycaseid+'_'+res
else:
    basecase=res

if (options.noad == False):
    cases.append(basecase+'_'+mymodel_adsp+'_ad_spinup')
if (options.nofn == False):
    cases.append(basecase+'_'+mymodel_fnsp)
if (options.notrans == False):
    cases.append(basecase+'_'+mymodel_trns)

for c in cases:
    model_startdate = 1
    if ('ad_spinup' in c):
        run_n_total = int(ny_ad)
    elif ('1850' in c):
        run_n_total = int(fsplen)
    elif ('20TR' in c or 'ICBCLM45' in c):
        run_n_total = int(translen)
        model_startdate = 1850
    else:
        run_n_total = int(fsplen)
    runblock =  min(int(options.runblock), run_n_total)
    n_submits = int(math.ceil(run_n_total / float(runblock)))


    mysubmit_type = 'qsub'
    if ('anvil' in options.machine or 'compy' in options.machine or 'cori' in options.machine):
        mysubmit_type = 'sbatch'
    #Create a .PBS site fullrun script to launch the full job 

    for n in range(0,n_submits):
        output = open('./temp/global_'+c+'_'+str(n)+'.pbs','w')
        if (os.path.isfile(caseroot+'/'+c+'/case.run')):
            input = open(caseroot+'/'+c+'/case.run')
        elif (os.path.isfile(caseroot+'/'+c+'/.case.run')):
            input = open(caseroot+'/'+c+'/.case.run')
        else:
            print 'case.run file not found.  Aborting'
            sys.exit(1)

        for s in input:
            if ("perl" in s or "python" in s):
                timestr=str(int(float(options.walltime)))+':'+str(int((float(options.walltime)- \
                            int(float(options.walltime)))*60))+':00'
                if (options.debug):
                    timestr='00:30:00'
                if ('cades' in options.machine):
                    output.write("#!/bin/bash -f\n")
                else:
                    output.write("#!/bin/csh -f\n")
                if (mysubmit_type == 'qsub'):
                    output.write('#PBS -l walltime='+timestr+'\n')
                else:
                    if (myproject != ''):
                      output.write('#SBATCH -A '+myproject+'\n')
                    output.write('#SBATCH --time='+timestr+'\n')
                    if ('anvil' in options.machine):
                      output.write('#SBATCH --partition=acme-centos6\n')
                      output.write('#SBATCH --account=condo\n')
                    if ('cori' in options.machine or 'edison' in options.machine):
                         if (options.debug):
                             output.write('#SBATCH --partition=debug\n')
                         else:
                             output.write('#SBATCH --partition=regular\n')
            elif ("#!" in s or "#PBS" in s or "#SBATCH" in s):
                output.write(s)
        input.close()
        output.write("\n")
   
        if (options.machine == 'cades'):
            output.write('source $MODULESHOME/init/bash\n')
            output.write('module unload python\n')
            output.write('module load python/2.7.12\n\n')
        if (options.machine == 'eos'):
            output.write('source $MODULESHOME/init/csh\n')
            output.write('module load nco\n')
            output.write('module unload python\n')
            output.write('module load python/2.7.5\n')
            output.write('module unload PrgEnv-intel\n')
            output.write('module load PrgEnv-gnu\n')
            output.write('module load python_numpy\n')
            output.write('module load python_scipy\n')
            output.write('module unload PrgEnv-gnu\n')
            output.write('module load PrgEnv-intel\n')
        if (options.machine == 'titan'):
            output.write('source $MODULESHOME/init/csh\n')
            output.write('module load nco\n')
            output.write('module load python\n')
            output.write('module load python_numpy/1.9.2\n')
            output.write('module load python_scipy/0.15.1\n')
        if ('cori' in options.machine or 'edison' in options.machine):
            output.write('source $MODULESHOME/init/csh\n')
            output.write('module unload python\n')
            output.write('module unload scipy\n')
            output.write('module unload numpy\n')
            output.write('module load python/2.7-anaconda\n')
            output.write('module load nco\n')                  
   
        this_run_n = runblock
        if (n == (n_submits-1)):
            this_run_n = run_n_total- (n * runblock)

        
        output.write("cd "+caseroot+'/'+c+"/\n")
        output.write("./xmlchange -id STOP_N -val "+str(this_run_n)+'\n')
        if (options.cplhist):
          output.write("./xmlchange -id REST_N -val 25\n")
        else:
          output.write("./xmlchange -id REST_N -val 20\n")   #Restart every 20 years in global mode
        output.write("./xmlchange -id RUN_STARTDATE -val "+(str(10000+model_startdate))[1:]+ \
                         '-01-01\n')                           
        if (n > 0):
            #change finidat 
            mylnd_in = open(caseroot+'/'+c+'/user_nl_clm','r')
            mylnd_out = open(caseroot+'/'+c+'/user_nl_clm_'+str(n), 'w')
            for s in mylnd_in:
                if ('finidat' in s):
                    mylnd_out.write(" finidat = '"+os.path.abspath(runroot)+"/"+c+"/run/"+c+ \
                              ".clm2.r."+str(10000+model_startdate)[1:]+"-01-01-00000.nc'\n")
                else:
                    mylnd_out.write(s)
            mylnd_in.close()
            mylnd_out.close()
            output.write(" cp user_nl_clm_"+str(n)+" user_nl_clm\n")
        model_startdate = model_startdate + runblock          
        output.write("./case.submit --no-batch\n")
        output.write("cd "+os.path.abspath(".")+'\n')
        if ('ad_spinup' in c and n == (n_submits-1)):
            if (options.bgc):
                output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+ \
                                 '/'+ad_case+'/run/ --casename '+ ad_case+' --restart_year '+ \
                             str(int(ny_ad)+1)+' --BGC\n')
            else:
                output.write("python adjust_restart.py --rundir "+os.path.abspath(runroot)+ \
                                 '/'+ad_case+'/run/ --casename '+ad_case+' --restart_year '+ \
                                 str(int(ny_ad)+1)+'\n')
        output.close()

        job_depend_run = submit('temp/global_'+c+'_'+str(n)+'.pbs',job_depend=job_depend_run, \
                                    submit_type=mysubmit_type)
        
