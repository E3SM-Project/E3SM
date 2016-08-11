#!/bin/bash
#
#
# Batch script to submit to create suite of ESMF mapping files
#
# Set up for yellowstone
# 
# yellowstone-specific batch commands:
#BSUB -P P00000000        # project number
#BSUB -n 16               # number of processors
#BSUB -R "span[ptile=16]" 
#BSUB -W 1:00             # wall-clock limit
#BSUB -q small            # queue
#BSUB -o regrid.%J.out    # ouput filename
#BSUB -e regrid.%J.err    # error filename
#BSUB -J gen_cesm_maps    # job name
#BSUB -N                  # send email upon job completion

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Set user-defined parameters here
#----------------------------------------------------------------------

griddir1="/glade/scratch/mlevy/grids"
griddir2="/glade/scratch/mlevy/grids"

fileocn="${griddir1}/tx0.1v2_090127.nc"
fileatm="${griddir1}/CESM/cseg/mapping/grids/fv0.9x1.25_070727.nc"
filertm="${griddir2}/SCRIPgrid_0.5x0.5_nomask_c110308.nc"
nameocn='tx0.1v2'
nameatm='fv0.9x1.25'
namertm='r0.5x0.5'

use_rtm=1
typeocn='global'
typeatm='global'

#----------------------------------------------------------------------
# Done setting user-defined parameters
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Stuff done in a machine-specific way
#----------------------------------------------------------------------

# Determine number of processors we're running on
host_array=($LSB_HOSTS)
REGRID_PROC=${#host_array[@]}

#----------------------------------------------------------------------
# Begin general script
#----------------------------------------------------------------------

cmdargs="--fileocn $fileocn --fileatm $fileatm --nameocn $nameocn --nameatm $nameatm --typeocn $typeocn --typeatm $typeatm"
if [ $use_rtm == 1 ]; then
  cmdargs="$cmdargs --filertm $filertm --namertm $namertm"
fi
cmdargs="$cmdargs --batch --nogridcheck"
env REGRID_PROC=$REGRID_PROC ./gen_cesm_maps.sh $cmdargs
