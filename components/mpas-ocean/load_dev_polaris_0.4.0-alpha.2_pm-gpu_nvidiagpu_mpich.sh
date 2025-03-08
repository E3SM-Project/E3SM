export POLARIS_BRANCH="/global/cfs/cdirs/e3sm/hgkang/polaris/polaris_MSdel2del4_nvidia"
export POLARIS_VERSION="0.4.0-alpha.2"

version_file="${POLARIS_BRANCH}/polaris/version.py"
code_version=$(cat $version_file)
if [[ "$code_version" != *"$POLARIS_VERSION"* ]]; then

echo "This load script is for a different version of polaris:"
echo "__version__ = '$POLARIS_VERSION'"
echo ""
echo "Your code is version:"
echo "$code_version"
echo ""
echo "You need to run ./configure_polaris_envs.py to update your conda "
echo "environment and load script."

else
# the right polaris version

echo Loading conda environment
source /global/cfs/projectdirs/e3sm/hgkang/miniconda3/etc/profile.d/conda.sh
conda activate dev_polaris_0.4.0-alpha.2
echo Done.
echo

if [[ -z "${NO_POLARIS_REINSTALL}" && -f "./setup.py" && \
-d "polaris" ]]; then
# safe to assume we're in the polaris repo
# update the polaris installation to point here
mkdir -p deploy_tmp/logs
echo Reinstalling polaris package in edit mode...
python -m pip install --no-deps --no-build-isolation -e . \
&> deploy_tmp/logs/install_polaris.log
echo Done.
echo
fi

echo Loading Spack environment...
source /global/cfs/cdirs/e3sm/software/polaris/pm-gpu/spack/dev_polaris_0_4_0_nvidiagpu_mpich/share/spack/setup-env.sh
spack env activate dev_polaris_0_4_0_nvidiagpu_mpich
module rm cray-hdf5-parallel &> /dev/null
module rm cray-netcdf-hdf5parallel &> /dev/null
module rm cray-parallel-netcdf &> /dev/null
module rm PrgEnv-gnu &> /dev/null
module rm PrgEnv-intel &> /dev/null
module rm PrgEnv-nvidia &> /dev/null
module rm PrgEnv-cray &> /dev/null
module rm PrgEnv-aocc &> /dev/null
module rm intel &> /dev/null
module rm intel-oneapi &> /dev/null
module rm cudatoolkit &> /dev/null
module rm climate-utils &> /dev/null
module rm matlab &> /dev/null
module rm craype-accel-nvidia80 &> /dev/null
module rm craype-accel-host &> /dev/null
module rm perftools-base &> /dev/null
module rm perftools &> /dev/null
module rm darshan &> /dev/null

module load PrgEnv-nvidia
module load nvidia/22.7
module load craype-x86-milan
module load libfabric/1.20.1
module load cudatoolkit/11.7
module load craype-accel-nvidia80
module load gcc-mixed/12.2.0
module load craype/2.7.20
module rm cray-mpich &> /dev/null
module load cray-mpich/8.1.25

module load cray-libsci/23.02.1.1

module rm cray-hdf5-parallel &> /dev/null
module rm cray-netcdf-hdf5parallel &> /dev/null
module rm cray-parallel-netcdf &> /dev/null
module load cray-hdf5-parallel/1.12.2.3
module load cray-netcdf-hdf5parallel/4.9.0.3
module load cray-parallel-netcdf/1.12.3.3

export NETCDF_PATH=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export NETCDF_C_PATH=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export NETCDF_FORTRAN_PATH=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export PNETCDF_PATH=$CRAY_PARALLEL_NETCDF_PREFIX

export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_MPIIO_DVS_MAXNODES=1
## purposefully omitting OMP variables that cause trouble in ESMF
# export OMP_STACKSIZE=128M
# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
export HDF5_USE_FILE_LOCKING=FALSE
## Not needed
# export PERL5LIB=/global/cfs/cdirs/e3sm/perl/lib/perl5-only-switch
export MPICH_GPU_SUPPORT_ENABLED=1

if [ -z "${NERSC_HOST:-}" ]; then
# happens when building spack environment
export NERSC_HOST="perlmutter"
fi
echo Done.
echo

export POLARIS_COMPILER=nvidiagpu
export POLARIS_MPI=mpich
export MPAS_EXTERNAL_LIBS=""
export NETCDF=${CRAY_NETCDF_HDF5PARALLEL_PREFIX}
export NETCDFF=${CRAY_NETCDF_HDF5PARALLEL_PREFIX}
export PNETCDF=${CRAY_PARALLEL_NETCDF_PREFIX}
#export PIO=/global/cfs/cdirs/e3sm/software/polaris/pm-gpu/spack/dev_polaris_0_4_0_nvidiagpu_mpich/var/spack/environments/dev_polaris_0_4_0_nvidiagpu_mpich/.spack-env/view
export PIO=/global/cfs/projectdirs/e3sm/hgkang/programs/ParallelIO_nvidia_new
export METIS_ROOT=/global/cfs/cdirs/e3sm/software/polaris/pm-gpu/spack/dev_polaris_0_4_0_nvidiagpu_mpich/var/spack/environments/dev_polaris_0_4_0_nvidiagpu_mpich/.spack-env/view
export PARMETIS_ROOT=/global/cfs/cdirs/e3sm/software/polaris/pm-gpu/spack/dev_polaris_0_4_0_nvidiagpu_mpich/var/spack/environments/dev_polaris_0_4_0_nvidiagpu_mpich/.spack-env/view
export PATH="/global/cfs/cdirs/e3sm/software/polaris/pm-gpu/spack/dev_polaris_soft_0_4_0/var/spack/environments/dev_polaris_soft_0_4_0/.spack-env/view/bin:$PATH"

export USE_PIO2=true
export OPENMP=true
export HDF5_USE_FILE_LOCKING=FALSE
export LOAD_POLARIS_ENV=/global/cfs/cdirs/e3sm/hgkang/polaris/polaris_MSdel2del4_nvidia/load_dev_polaris_0.4.0-alpha.2_pm-gpu_nvidiagpu_mpich.sh
export POLARIS_MACHINE=pm-gpu

# the right polaris version
fi
