#!/bin/bash
#
#
# Batch script to submit to create ESMF mapping file
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
#BSUB -J create_ESMF_map  # job name
#BSUB -N                  # send email upon job completion

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Set user-defined parameters here
#----------------------------------------------------------------------

griddir="/glade/scratch/mlevy/grids"

filesrc="$griddir/tx0.1v2_090127.nc"
filedst="$griddir/fv0.9x1.25_070727.nc"
namesrc='tx0.1v2'
namedst='fv0.9x1.25'

typesrc='global'
typedst='global'
maptype='aave'

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

cmdargs="--filesrc $filesrc --filedst $filedst --namesrc $namesrc --namedst $namedst --typesrc $typesrc --typedst $typedst --maptype $maptype --batch"
env REGRID_PROC=$REGRID_PROC ./create_ESMF_map.sh $cmdargs
