#!/usr/bin/env python
# script for running E3SM-MMF benchmark simulations using the 2020 INICTE allocation (CLI115)
# Branch for this campaign: https://github.com/E3SM-Project/E3SM/tree/whannah/incite-2020
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import os, numpy as np, subprocess as sp, datetime
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

acct = 'cli115'
case_dir = os.getenv('HOME')+'/E3SM/Cases/'
src_dir  = os.getenv('HOME')+'/E3SM/E3SM_SRC1/'

### flags to control config/build/submit sections
# clean        = True
newcase      = True
config       = True
build        = True
submit       = True
# continue_run = True

### run duration and resubmission (remember to set continue_run)
stop_opt,stop_n,resub = 'ndays',1,0

### same settings for all runs
compset = 'F-MMFXX'
ne,npg  = 45,2
grid    = f'ne{ne}pg{npg}_r05_oECv3'
# grid    = f'ne{ne}pg{npg}_ne{ne}pg{npg}'
crm_dx  = 3200 
crm_dt  = 10
rad_nx  = 2

### time stamps for distinguishing different settings
# timestamp = '20201012' # same config as larger INCITE runs, task_per_node=84
timestamp = '20201013' # change to 1 day run, task_per_node=84
# timestamp = '20201014' # change CPU task_per_node=48 (also fix how other component tasks are set)
# timestamp = '20201015' # change CPU task_per_node=12

### case specifics
# arch='GNUCPU'
# arch='GNUGPU'

### CPU
arch='GNUCPU'
num_nodes = 578
# num_nodes = 290
# num_nodes = 145
# num_nodes = 73

##GPU
# arch='GNUGPU'
# num_nodes = 1013
# num_nodes = 507
# num_nodes = 254
# num_nodes = 127
# num_nodes = 64


# num_nodes = 32
# num_nodes = 16
# num_nodes = 8
# num_nodes = 4
# num_nodes = 2
# num_nodes = 1

### specify full case name
case_list_common = ['INCITE2020','BENCHMARK',arch,grid,compset,'NLEV_50','CRMNX_32']
# case_list = case_list_common+['CRMNY_32','MOMFB']
case_list = case_list_common+['CRMNY_32']
# case_list = case_list_common # 2D

case='.'.join(case_list+[f'NN_{num_nodes}',timestamp] )

# Impose wall limits for Summits
# if num_nodes>=  1: walltime =  '2:00'
# if num_nodes>= 46: walltime =  '6:00'
# if num_nodes>= 92: walltime = '12:00'
# if num_nodes>=922: walltime = '24:00'
walltime = '6:00'

### add these modifiers to enable debug mode or state variable checks
# case += '.debug-on'
# case += '.checks-on'

### specify atmos initial condition file
# init_file_dir = '/gpfs/alpine/scratch/hannah6/cli115/e3sm_scratch/init_files'
init_file_dir = '/gpfs/alpine/scratch/hannah6/cli115/HICCUP/data/'
params = [p.split('_') for p in case.split('.')]
for p in params:
   if p[0]=='NLEV': 
      if p[1]!='72': init_file_atm = f'HICCUP.cami_mam3_Linoz_ne{ne}np4.L{p[1]}.nc'

### specify land initial condition file
if ne==45:
   land_init_path = '/gpfs/alpine/scratch/hannah6/cli115/e3sm_scratch/init_files'
   land_init_file = 'CLM_spinup.ICRUELM.ne45pg2_r05_oECv3.20-yr.2010-10-01.elm.r.2006-01-01-00000.nc'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

dtime = 20*60   # GCM physics time step

if 'dtime' in locals(): ncpl  = 86400 / dtime

num_dyn = ne*ne*6

if 'task_per_node' not in locals():
   if arch=='GNUCPU': task_per_node = 84
   if arch=='GNUGPU': task_per_node = 12
#-------------------------------------------------------------------------------
# Define run command
#-------------------------------------------------------------------------------
# Set up terminal colors
class tcolor:
   ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd,suppress_output=False,execute=True):
   if suppress_output : cmd = cmd + ' > /dev/null'
   msg = tcolor.GREEN + cmd + tcolor.ENDC
   print(f'\n{msg}')
   if execute: os.system(cmd)
   return
#---------------------------------------------------------------------------------------------------
# Create new case
#---------------------------------------------------------------------------------------------------
if newcase :
   cmd = src_dir+'cime/scripts/create_newcase --case '+case_dir+case
   cmd = cmd + ' --compset '+compset+' --res '+grid
   # cmd = cmd + ' --pecount '+str(num_nodes*task_per_node)+'x1 '
   if arch=='GNUCPU': cmd = cmd + ' -compiler gnu    '
   if arch=='GNUGPU': cmd = cmd + ' -compiler gnugpu '
   run_cmd(cmd)

   # Change run directory to be next to bld directory
   os.chdir(case_dir+case+'/')
   memberwork = os.getenv('MEMBERWORK')
   run_cmd(f'./xmlchange -file env_run.xml RUNDIR=\'{memberwork}/{acct}/e3sm_scratch/{case}/run\' ' )
   
   run_cmd(f'./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE    -val {task_per_node} ')
   run_cmd(f'./xmlchange -file env_mach_pes.xml -id MAX_MPITASKS_PER_NODE -val {task_per_node} ')

#---------------------------------------------------------------------------------------------------
# Configure
#---------------------------------------------------------------------------------------------------
os.chdir(case_dir+case+'/')
if config : 
   #-------------------------------------------------------
   if arch=='GNUGPU': 
      pcols = np.ceil( (ne**2*6*npg**2) / (num_nodes*task_per_node) )
      run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -pcols {int(pcols)} \" ' )
   #-------------------------------------------------------
   if 'init_file_atm' in locals():
      file = open('user_nl_eam','w')
      file.write(f' ncdata = \'{init_file_dir}/{init_file_atm}\'\n')
      file.close()
   #-------------------------------------------------------
   # Specify CRM and RAD columns
   params = [p.split('_') for p in case.split('.')]
   for p in params:
      if p[0]=='CRMNX': run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -crm_nx {p[1]} \" ')
      if p[0]=='CRMNY': run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -crm_ny {p[1]} \" ')
      if p[0]=='NLEV' and p[1] != '72' : 
         nlev = p[1]; crm_nz = None
         if nlev== '50': crm_nz =  '46'
         if nlev=='100': crm_nz =  '95'
         if nlev=='120': crm_nz = '115'
         if crm_nz is None: raise ValueError('No value of crm_nz specified')
         run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -nlev {p[1]} -crm_nz {crm_nz} \" ')
   #-------------------------------------------------------
   # add in the settings that are specified above
   rad_ny = rad_nx if 'CRMNY' in case else 1
   run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -crm_dx {crm_dx} -crm_dt {crm_dt}  \" ')
   run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -crm_nx_rad {rad_nx} -crm_ny_rad {rad_ny} \" ')
   #-------------------------------------------------------
   # Add special MMF options based on case name
   cpp_opt = ''
   if '.MOMFB.'  in case: cpp_opt += ' -DMMF_MOMENTUM_FEEDBACK'
   if 'debug-on' in case: cpp_opt += ' -DYAKL_DEBUG'

   if cpp_opt != '' :
      cmd  = f'./xmlchange --append -file env_build.xml -id CAM_CONFIG_OPTS'
      cmd += f' -val \" -cppdefs \'-DMMF_DIR_NS {cpp_opt} \'  \" '
      run_cmd(cmd)
   #-------------------------------------------------------
   # Set tasks and threads
   atm_ntasks = task_per_node*num_nodes
   cmd = './xmlchange -file env_mach_pes.xml '
   cmd += f' NTASKS_ATM={atm_ntasks}'
   alt_ntask = task_per_node * min([ne,num_nodes])
   cmd += f',NTASKS_LND={alt_ntask},NTASKS_CPL={alt_ntask}'
   cmd += f',NTASKS_OCN={alt_ntask},NTASKS_ICE={alt_ntask}'
   cmd += f',NTASKS_ROF=1,NTASKS_WAV=1,NTASKS_GLC=1,NTASKS_ESP=1,NTASKS_IAC=1'
   run_cmd(cmd)
   # run_cmd('./xmlchange -file env_mach_pes.xml NTHRDS_ATM=2,NTHRDS_CPL=2,NTHRDS_LND=1')
   #-------------------------------------------------------
   # 64_data format is needed for large numbers of columns (GCM or CRM)
   run_cmd('./xmlchange PIO_NETCDF_FORMAT=\"64bit_data\" ')
   #-------------------------------------------------------
   # Run case setup
   if clean : run_cmd('./case.setup --clean')
   run_cmd('./case.setup --reset')
#---------------------------------------------------------------------------------------------------
# Build
#---------------------------------------------------------------------------------------------------
if build : 
   # run_cmd('./xmlchange PIO_VERSION=1')
   if 'debug-on' in case : run_cmd('./xmlchange -file env_build.xml -id DEBUG -val TRUE ')
   if clean : run_cmd('./case.build --clean')
   run_cmd('./case.build')
#---------------------------------------------------------------------------------------------------
# Write the namelist options and submit the run
#---------------------------------------------------------------------------------------------------
if submit : 
   # Change inputdata from default due to permissions issue
   run_cmd('./xmlchange DIN_LOC_ROOT=/gpfs/alpine/cli115/scratch/hannah6/inputdata ')
   
   #-------------------------------------------------------
   # First query some stuff about the case
   #-------------------------------------------------------
   (din_loc_root, err) = sp.Popen('./xmlquery DIN_LOC_ROOT    -value', \
                                     stdout=sp.PIPE, shell=True, universal_newlines=True).communicate()
   (config_opts, err) = sp.Popen('./xmlquery CAM_CONFIG_OPTS -value', \
                                     stdout=sp.PIPE, shell=True, universal_newlines=True).communicate()
   config_opts = ' '.join(config_opts.split())   # remove extra spaces to simplify string query
   ntasks_atm = None
   (ntasks_atm     , err) = sp.Popen('./xmlquery NTASKS_ATM    -value', \
                                     stdout=sp.PIPE, shell=True, universal_newlines=True).communicate()
   ntasks_atm = float(ntasks_atm)
   #-------------------------------------------------------
   # Namelist options
   #-------------------------------------------------------
   nfile = 'user_nl_eam'
   file = open(nfile,'w') 
   #------------------------------
   # Other namelist stuff
   #------------------------------
   # file.write(' srf_flux_avg = 1 \n')              # sfc flux smoothing (for MMF stability)
   # file.write(f' crm_accel_factor = 3 \n')         # CRM acceleration factor (default is 2)

   if num_dyn<ntasks_atm: file.write(' dyn_npes = '+str(num_dyn)+' \n')   # limit dynamics tasks

   if 'checks-on' in case: file.write(' state_debug_checks = .true. \n')

   # disable cami file output
   file.write(" inithist = \'NONE\' \n") # ENDOFRUN / NONE

   if 'init_file_atm' in locals():
      file.write(f' ncdata = \'{init_file_dir}/{init_file_atm}\'\n')

   file.close()
   #-------------------------------------------------------
   # CLM namelist
   #-------------------------------------------------------
   if 'land_init_file' in locals():
      nfile = 'user_nl_elm'
      file = open(nfile,'w')
      file.write(f' finidat = \'{land_init_path}/{land_init_file}\' \n')
      file.close()
   #-------------------------------------------------------
   # Set some run-time stuff
   #-------------------------------------------------------
   if 'ncpl' in locals(): run_cmd(f'./xmlchange ATM_NCPL={str(ncpl)}')
   run_cmd(f'./xmlchange STOP_OPTION={stop_opt},STOP_N={stop_n},RESUBMIT={resub}')
   run_cmd(f'./xmlchange JOB_QUEUE=batch,JOB_WALLCLOCK_TIME={walltime}')
   run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')

   # disable restart files
   run_cmd(f'./xmlchange REST_OPTION=never')

   # # Restart Frequency
   # run_cmd('./xmlchange -file env_run.xml REST_OPTION={stop_opt},REST_N={stop_n}')

   if continue_run :
      run_cmd('./xmlchange -file env_run.xml CONTINUE_RUN=TRUE ')   
   else:
      run_cmd('./xmlchange -file env_run.xml CONTINUE_RUN=FALSE ')
   #-------------------------------------------------------
   # Submit the run
   #-------------------------------------------------------
   run_cmd('./case.submit')
#---------------------------------------------------------------------------------------------------
# Print the case name again
#---------------------------------------------------------------------------------------------------
print(f'\n  case : {case}\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
