export POLARIS_BRANCH="/autofs/nccs-svm1_home1/mpetersen/repos/polaris/main"
export POLARIS_VERSION="0.5.0-alpha.3"

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
source /ccs/home/mpetersen/miniforge3/etc/profile.d/conda.sh
conda activate dev_polaris_0.5.0-alpha.3
echo Done.
echo

if [[ -z "${NO_POLARIS_REINSTALL}" && -f "./pyproject.toml" && \
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
source /ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_5_0_crayclanggpu_mpich/share/spack/setup-env.sh
spack env activate dev_polaris_0_5_0_crayclanggpu_mpich
module reset >& /dev/null
module switch Core/24.07 >& /dev/null
module switch PrgEnv-cray PrgEnv-cray/8.3.3 >& /dev/null
module switch cce cce/18.0.1 >& /dev/null
module switch craype craype/2.7.20 >& /dev/null
module load craype-accel-amd-gfx90a
module load rocm/6.2.4

module load cray-libsci/22.12.1.1

module load cray-hdf5-parallel/1.12.2.11
module load cray-netcdf-hdf5parallel/4.9.0.11
module load cray-parallel-netcdf/1.12.3.11

export NETCDF_C_PATH=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export NETCDF_FORTRAN_PATH=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export PNETCDF_PATH=$CRAY_PARALLEL_NETCDF_PREFIX

export HDF5_USE_FILE_LOCKING=FALSE
echo Done.
echo

export POLARIS_COMPILER=crayclanggpu
export POLARIS_MPI=mpich
export MPAS_EXTERNAL_LIBS=""
export NETCDF=$(dirname $(dirname $(which nc-config)))
export NETCDFF=$(dirname $(dirname $(which nf-config)))
export PNETCDF=$(dirname $(dirname $(which pnetcdf-config)))
export PIO=/ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_5_0_crayclanggpu_mpich/var/spack/environments/dev_polaris_0_5_0_crayclanggpu_mpich/.spack-env/view
export METIS_ROOT=/ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_5_0_crayclanggpu_mpich/var/spack/environments/dev_polaris_0_5_0_crayclanggpu_mpich/.spack-env/view
export PARMETIS_ROOT=/ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_5_0_crayclanggpu_mpich/var/spack/environments/dev_polaris_0_5_0_crayclanggpu_mpich/.spack-env/view
export PATH="/ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_soft_0_5_0/var/spack/environments/dev_polaris_soft_0_5_0/.spack-env/view/bin:$PATH"

export USE_PIO2=true
export OPENMP=false
export HDF5_USE_FILE_LOCKING=FALSE
export LOAD_POLARIS_ENV=/autofs/nccs-svm1_home1/mpetersen/repos/polaris/main/load_dev_polaris_0.5.0-alpha.2_frontier_crayclanggpu_mpich.sh
export POLARIS_MACHINE=frontier

# the right polaris version
fi
