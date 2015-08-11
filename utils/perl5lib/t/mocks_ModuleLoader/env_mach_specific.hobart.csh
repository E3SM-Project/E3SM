#!/usr/bin/env csh -f 
#===============================================================================
# Automatically generated module settings for hobart
#===============================================================================

source /usr/share/Modules/init/csh
module purge 
if ( $COMPILER == "intel" ) then
	module load compiler/intel/15.0.2.164
endif
if ( $COMPILER == "intel" && $MPILIB == "mvapich2" ) then
	module unload mpi/intel/openmpi-1.8.1-qlc
	module load mpi/intel/mvapich2-1.8.1-qlc
endif
if ( $COMPILER == "pgi" ) then
	module load compiler/pgi/15.1
endif
if ( $COMPILER == "pgi" && $MPILIB == "mvapich2" ) then
	module unload mpi/pgi/openmpi-1.8.1-qlc
	module load mpi/pgi/mvapich2-1.8.1-qlc
endif
if ( $COMPILER == "nag" ) then
	module load compiler/nag/6.0
	./xmlchange MPILIB=openmpi
endif
if ( $COMPILER == "gnu" ) then
	module load compiler/gnu/4.8.3
endif
