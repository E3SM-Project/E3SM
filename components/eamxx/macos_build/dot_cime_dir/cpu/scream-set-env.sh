#! /bin/bash

export SCREAM_INPUT_DIR="${HOME}/.cime/scream-input"
export SCREAM_BASELINES_DIR="${HOME}/.cime/scream-data/master-baselines"

# cee_modpath="/projects/aue/cee/deploy/9d138342/linux-rhel8-x86_64/gcc-12.3.0"
export NetCDF_C_ROOT="/opt/homebrew/Cellar/netcdf/4.9.3"
export PnetCDF_C_ROOT="/opt/homebrew/Cellar/pnetcdf/1.14.0"
export NetCDF_Fortran_ROOT="/opt/homebrew/Cellar/netcdf-fortran/4.6.2"
export HDF5_ROOT="/opt/homebrew/Cellar/hdf5/1.14.6"
# export LAPACK_ROOT="${cee_modpath}/netlib-lapack-3.11.0-xyhaan3"
# export MPI_ROOT="${cee_modpath}/openmpi-4.1.6-hfp5vwq"

export LD_LIBRARY_PATH="${NetCDF_C_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PnetCDF_C_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${NetCDF_Fortran_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}"
# export LD_LIBRARY_PATH="${LAPACK_ROOT}/lib64:${LD_LIBRARY_PATH}"

export LIBRARY_PATH="${NetCDF_C_ROOT}/lib:${LIBRARY_PATH}"
export LIBRARY_PATH="${PnetCDF_C_ROOT}/lib:${LIBRARY_PATH}"
export LIBRARY_PATH="${NetCDF_Fortran_ROOT}/lib:${LIBRARY_PATH}"
export LIBRARY_PATH="${HDF5_ROOT}/lib:${LIBRARY_PATH}"
# export LIBRARY_PATH="${LAPACK_ROOT}/lib64:${LIBRARY_PATH}"

export PKG_CONFIG_PATH="${NetCDF_C_ROOT}/lib/pkgconfig:${PKG_CONFIG_PATH}"
export PKG_CONFIG_PATH="${PnetCDF_C_ROOT}/lib/pkgconfig:${PKG_CONFIG_PATH}"
export PKG_CONFIG_PATH="${NetCDF_Fortran_ROOT}/lib/pkgconfig:${PKG_CONFIG_PATH}"
export PKG_CONFIG_PATH="${HDF5_ROOT}/lib/pkgconfig:${PKG_CONFIG_PATH}"
# export PKG_CONFIG_PATH="${LAPACK_ROOT}/lib64/pkgconfig:${PKG_CONFIG_PATH}"

export PATH="${NetCDF_C_ROOT}/bin:${PATH}"
export PATH="${PnetCDF_C_ROOT}/bin:${PATH}"
export PATH="${NetCDF_Fortran_ROOT}/bin:${PATH}"
export PATH="${HDF5_ROOT}/bin:${PATH}"
# export PATH="${MPI_ROOT}/bin:${PATH}"
