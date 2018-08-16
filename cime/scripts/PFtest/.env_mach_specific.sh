module purge 
module load PE-gnu mkl/2017 cmake/3.6.1 python/2.7.12 nco/4.6.4 hdf5-parallel/1.8.17 netcdf-hdf5parallel/4.3.3.1
export PETSC_PATH=/software/user_tools/current/cades-ccsi/petsc4pf/openmpi-1.10-gcc-5.3
export CLM_PFLOTRAN_COUPLED=FALSE
export CLM_PFLOTRAN_COLMODE=FALSE
export CLM_PFLOTRAN_SOURCE_DIR=/lustre/or-hydra/cades-ccsi/omearata/models/pflotran-interface/src/clm-pflotran