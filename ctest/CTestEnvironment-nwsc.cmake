#==============================================================================
#
#  This file sets the environment variables needed to configure and build
#  CTest on the NCAR Wyoming Supercomputing Center (yellowstone).
#
#==============================================================================

# Load modules
execute_process (COMMAND module save
                 COMMAND module reset
                 COMMAND module swap intel intel/15.0.3
                 COMMAND module load git/2.3.0
                 COMMAND module load cmake/3.0.2
                 COMMAND module load netcdf-mpi/4.3.3.1
                 COMMAND module load pnetcdf/1.6.0)
                 
# Assume all package locations (NetCDF, PnetCDF, etc) are already
# set with existing environment variables: NETCDF, PNETCDF, etc.

set (ENV{CC} "mpicc")
set (ENV{FC} "mpif90")
