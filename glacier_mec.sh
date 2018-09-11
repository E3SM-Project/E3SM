#!/bin/bash

#######################################################################
#######################################################################
#######  Script to run ACME in ensemble mode using CIME's multi instance feature
######   Currently runs a test FC5 with 2 active cases
#######################################################
#######  BEGIN USER DEFINED SETTINGS

# Set the name of your case here
casename=ICLM45GLCMEC

# Set the case directory here
casedirectory=$PROJWORK/cli115/$USER
 
# Name of machine you are running on (i.e. edison, anvil, etc)                                                    
machine=titan

# Name of project to run on, if submitting to queue
projectname=cli115

# Name of run type
# compset=ICLM45GLCMEC         # 2000_DATM%QIA_CLM45%CN_SICE_SOCN_MOSART_CISM1_SWAV_TEST
# compset=IG1850CLM45          # 1850_DATM%QIA_CLM45%SP_SICE_SOCN_MOSART_CISM1_SWAV
# compset=IGCLM45              # 2000_DATM%QIA_CLM45%SP_SICE_SOCN_MOSART_CISM1_SWAV
# compset=IG20TRCLM45          # 20TR_DATM%QIA_CLM45%SP_SICE_SOCN_MOSART_CISM1_SWAV
# compset=IGRCP85CLM45CN       # RCP8_DATM%QIA_CLM45%CN_SICE_SOCN_MOSART_CISM1_SWAV
# compset=IGRCP45CLM45CN       # RCP4_DATM%QIA_CLM45%CN_SICE_SOCN_MOSART_CISM1_SWAV

# compset=IGCLM45IS2           # 2000_DATM%QIA_CLM45%SP_SICE_SOCN_MOSART_CISM2P_SWAV



compset=IGCLM45_MLI          # 2000_DATM%QIA_CLM45%SP_SICE_SOCN_MOSART_MALI%SIA_SWAV
# compset=TG_MLI               # 2000_SATM_DLND%SCPL_SICE_SOCN_SROF_MALI_SWAV
# compset=MALISIA              # RCP8_SATM_SLND_SICE_SOCN_SROF_MALI%SIA_SWAV
# compset=MALI                 # RCP8_SATM_SLND_SICE_SOCN_SROF_MALI_SWAV
# grid=f09_g16_g
grid=f09_g16_a

# compset=A_CRYO               # 2000_CAM5%AV1C-L_CLM45%SPBC_MPASCICE_MPASO_MOSART_MALI_SWAV
# compset=A_BG1850CN           # 1850_CAM5_CLM45%CN_MPASCICE_MPASO_MOSART_MALI%SIA_SWAV
# compset=MPAS_LISIO_TEST      # 2000_DATM%NYF_SLND_MPASCICE_MPASO_DROF%NYF_MALI%SIA_SWAV
# grid=f09_g16_g
# grid=f09_g16_a
# grid=T62_m120_g
# grid=f09_oEC_a
# grid=T62_oQU120_ais20


# grid resolution
# grid=f09_g16_g
#      alias: f09_g16_g (only for compsets that are _MALI )
#        non-default grids are: atm:0.9x1.25  lnd:0.9x1.25  ocnice:gx1v6  rof:r05  glc:mpas.gis20km  wav:null
#        mask is: gx1v6
# grid=f09_g16_a
#      alias: f09_g16_a (only for compsets that are _MALI )
#        non-default grids are: atm:0.9x1.25  lnd:0.9x1.25  ocnice:gx1v6  rof:r05  glc:mpas.ais20km  wav:null
#        mask is: gx1v6
#
# grid=f09_oEC_a
#      alias: f09_oEC_a (only for compsets that are MPASO.*_MALI )
#        non-default grids are: atm:0.9x1.25  lnd:0.9x1.25  ocnice:oEC60to30  rof:r05  glc:mpas.ais20km  wav:null
#        mask is: oEC60to30
# grid=T62_m120_g
#      alias: T62_m120_g (only for compsets that are MPASO.*_MALI )
#        non-default grids are: atm:T62  lnd:T62  ocnice:mpas120  rof:rx1  glc:mpas.gis20km  wav:null
#        mask is: mpas120
# grid=T62_oQU120_ais20
#      alias: T62_oQU120_ais20 (only for compsets that are MPASO.*_MALI )
#        non-default grids are: atm:T62  lnd:T62  ocnice:oQU120  rof:rx1  glc:mpas.ais20km  wav:null
#        mask is: oQU120


# User enter any needed modules to load or use below
module load python/2.7.9

# Directory where code lives
code_dir=$HOME/ACME

# Code tag name
code_tag=ACME

# Create new case
$code_dir/$code_tag/cime/scripts/create_newcase -case $casedirectory/$casename -mach $machine -project $projectname -compset $compset -res $grid

cd $casedirectory/$casename

./pelayout 
 
./xmlquery TOTALPES

# Enter CAM namelist options 
#cat <<EOF >> user_nl_cam
# fincl1 = 'SWCF:A','LWCF:A'
# hist_nhtfrq = 0,-24
#EOF
# Enter CLM namelist options 
cat <<EOF >> user_nl_clm
 flndtopo = '/lustre/atlas1/cli900/world-shared/cesm/inputdata/lnd/clm2/griddata/topodata_0.9x1.25_USGS_070110.nc'
 fglcmask = '/lustre/atlas1/cli900/world-shared/cesm/inputdata/lnd/clm2/griddata/glcmaskdata_0.9x1.25_GIS_AIS.nc'
EOF
cat <<EOF >> user_nl_mali
 mali_use_albany='FALSE'
EOF

# set up the case
  ./case.setup

# Build the case 
  ./case.build

  ./xmlchange STOP_N=32
  ./xmlchange JOB_WALLCLOCK_TIME=02:30:00

# Submit case 
  ./case.submit
  
