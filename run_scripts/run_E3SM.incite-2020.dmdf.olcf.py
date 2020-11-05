#!/usr/bin/env python
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

# same settings for all runs
arch,compset = 'GNUCPU','F-MMF1-AQP1'
ne,npg  = 4,2
grid    = f'ne{ne}pg{npg}_ne{ne}pg{npg}'

num_nodes=4
if 'GPU' in arch : task_per_node=12
if 'CPU' in arch : task_per_node=48

crm_nx, crm_ny, crm_dx, crm_dt, rad_nx, nlev, crm_nz = 32, 1, 3200, 10, 2, 50, 46  # 2D
# crm_nx, crm_ny, crm_dx, crm_dt, rad_nx, nlev, crm_nz = 32, 32, 3200, 10, 2, 50, 46  # 3D

timestamp = '20201104' # new batch with timestep output

# common parts of the case name
case_list = ['INCITE2020','DMDF',arch,grid,compset,f'NLEV_{nlev}',f'CRMNX_{crm_nx}']

# Impose wall limits for Summits
if num_nodes>=  1: walltime =  '2:00'
if num_nodes>= 46: walltime =  '6:00'
if num_nodes>= 92: walltime = '12:00'
if num_nodes>=922: walltime = '24:00'
# walltime = '2:00'

# specify atmos initial condition file
init_file_dir = '/gpfs/alpine/scratch/hannah6/cli115/HICCUP/data/'
init_file_atm = f'HICCUP.AQUA.ne{ne}np4.L{nlev}.nc'
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

num_dyn = ne*ne*6

dtime = 20*60   # GCM physics time step
if 'dtime' in locals(): ncpl  = 86400 / dtime

if 'task_per_node' not in locals():
   if arch=='GNUCPU': task_per_node = 84
   if arch=='GNUGPU': task_per_node = 12

atm_ntasks = task_per_node*num_nodes
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

#-------------------------------------------------------------------------------
# Create new case
#-------------------------------------------------------------------------------
if newcase :
   cmd = f'{src_dir}/cime/scripts/create_newcase --case {case_dir}{case}'
   cmd += f' --compset {compset} --res {grid} --pecount {atm_ntasks}x1 '
   if arch=='GNUCPU': cmd = cmd + ' -compiler gnu    '
   if arch=='GNUGPU': cmd = cmd + ' -compiler gnugpu '
   run_cmd(cmd)

   # Change run directory to be next to bld directory
   os.chdir(case_dir+case+'/')

   run_cmd(f'./xmlchange -file env_run.xml RUNDIR=\''+os.getenv('MEMBERWORK')+f'/{acct}/e3sm_scratch/{case}/run\' ' )
   run_cmd(f'./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE    -val {task_per_node} ')
   run_cmd(f'./xmlchange -file env_mach_pes.xml -id MAX_MPITASKS_PER_NODE -val {task_per_node} ')

#-------------------------------------------------------------------------------
# Configure
#-------------------------------------------------------------------------------
os.chdir(case_dir+case+'/')
if config : 
   
   # adjust pcols for GPU cases
   if arch=='GNUGPU': 
      pcols = np.ceil( (ne**2*6*npg**2) / (num_nodes*task_per_node) )
      run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -pcols {int(pcols)} \" ' )
   
   # if using non-standard atmos IC set it here to avoid config error
   if 'init_file_atm' in locals():
      file = open('user_nl_eam','w')
      file.write(f' ncdata = \'{init_file_dir}/{init_file_atm}\'\n')
      file.close()

   # Specify CRM and RAD settings
   rad_ny = rad_nx if crm_ny==1 else 1
   run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -crm_dx {crm_dx} -crm_dt {crm_dt}  \" ')
   run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -crm_nx_rad {rad_nx} -crm_ny_rad {rad_ny} \" ')
   run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -crm_nx {crm_nx} -crm_ny {crm_ny} \" ')
   run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -nlev {nlev} -crm_nz {crm_nz} \" ')
   
   # Set tasks and threads
   if arch == 'GNUGPU':
      cmd = './xmlchange -file env_mach_pes.xml '
      cmd += f' NTASKS_ATM={atm_ntasks}'
      if ne==4:  alt_ntask = task_per_node*4
      cmd += f',NTASKS_LND={alt_ntask},NTASKS_CPL={alt_ntask}'
      cmd += f',NTASKS_OCN={alt_ntask},NTASKS_ICE={alt_ntask}'
      cmd += f',NTASKS_ROF={task_per_node},NTASKS_WAV={task_per_node},NTASKS_GLC={task_per_node}'
      cmd += f',NTASKS_ESP=1,NTASKS_IAC=1'
      run_cmd(cmd)
      # run_cmd('./xmlchange -file env_mach_pes.xml NTHRDS_ATM=2,NTHRDS_CPL=2,NTHRDS_LND=1')
   
   # 64_data format is needed for large numbers of columns (GCM or CRM)
   run_cmd('./xmlchange PIO_NETCDF_FORMAT=\"64bit_data\" ')
   
   # Run case setup
   if clean : run_cmd('./case.setup --clean')
   run_cmd('./case.setup --reset')

#-------------------------------------------------------------------------------
# Build
#-------------------------------------------------------------------------------
if build : 
   if clean : run_cmd('./case.build --clean')
   run_cmd('./case.build')

#-------------------------------------------------------------------------------
# Write the namelist options and submit the run
#-------------------------------------------------------------------------------
if submit : 
   # Change inputdata from default due to permissions issue
   run_cmd('./xmlchange DIN_LOC_ROOT=/gpfs/alpine/cli115/scratch/hannah6/inputdata ')
   
   #-------------------------------------------------------
   # Namelist options
   nfile = 'user_nl_eam'
   file = open(nfile,'w') 
   atm_ntasks = task_per_node*num_nodes
   if num_dyn<atm_ntasks: file.write(' dyn_npes = '+str(num_dyn)+' \n')   # limit dynamics tasks
   if 'init_file_atm' in locals(): file.write(f' ncdata = \'{init_file_dir}/{init_file_atm}\'\n')
   file.close()

   #-------------------------------------------------------
   # Set some run-time stuff and submit
   if 'ncpl' in locals(): run_cmd(f'./xmlchange ATM_NCPL={str(ncpl)}')
   run_cmd(f'./xmlchange STOP_OPTION={stop_opt},STOP_N={stop_n},RESUBMIT={resub}')
   run_cmd(f'./xmlchange JOB_QUEUE=batch,JOB_WALLCLOCK_TIME={walltime}')
   run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')

   if continue_run :
      run_cmd('./xmlchange -file env_run.xml CONTINUE_RUN=TRUE ')   
   else:
      run_cmd('./xmlchange -file env_run.xml CONTINUE_RUN=FALSE ')
   
   run_cmd('./case.submit')

#-------------------------------------------------------------------------------
# Print the case name again
#-------------------------------------------------------------------------------
print(f'\n  case : {case}\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
