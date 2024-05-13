from utils import expect, run_cmd_no_fail, ensure_psutil

import os, sys, pathlib
ensure_psutil()
import psutil

CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","cime")

# MACHINE -> (env_setup,                      # list of shell commands to set up scream-approved env
#             compilers,                      # list of compilers [CXX, F90, C]
#             batch submit prefix,            # string shell commmand prefix
#             pre-existing baselines root dir # string path)
MACHINE_METADATA = {
    # NOTE: blake must use a different minor version of python than weaver
    # this is so pip installs do not go into the same location for both machines.
    # We have found some python modules are not shareable between these two machinse.
    "blake"    : (["module purge", "module load python/3.8.8/gcc/10.2.0", "module unload gcc/10.2.0", "module load openmpi/2.1.2 zlib git/2.9.4 cmake/3.19.3",
                   "export PATH=/ascldap/users/projects/e3sm/scream/libs/netcdf-fortran/install/blake/bin:$PATH",
                   "export PATH=/ascldap/users/projects/e3sm/scream/libs/netcdf-c/install/blake/bin:$PATH",
                   "export PATH=/ascldap/users/projects/e3sm/scream/libs/pnetcdf/install/blake/bin:$PATH",
                   "export INTEL_LICENSE_FILE=/home/projects/x86-64/intel/licenses/USE_SERVER.lic"  # speeds up license lookup
                  ],
                 ["mpicxx","mpifort","mpicc"],
                  "salloc -N 1 srun -n1 --preserve-env",
                  "/home/projects/e3sm/scream/pr-autotester/master-baselines/blake/"),
    "weaver"   : (["source /etc/profile.d/modules.sh", "module purge", "module load cmake/3.25.1 git/2.39.1 python/3.10.8 py-netcdf4/1.5.8 gcc/11.3.0 cuda/11.8.0 openmpi netcdf-c netcdf-fortran parallel-netcdf netlib-lapack",
                 ],
                 ["mpicxx","mpifort","mpicc"],
                  "bsub -I -q rhel8 -n 4 -gpu num=4",
                  "/home/projects/e3sm/scream/pr-autotester/master-baselines/weaver/"),
    "mappy"   : (["module purge", "module load sems-archive-env acme-env acme-cmake/3.26.3 acme-gcc/11.2.0 sems-archive-git/2.10.1 acme-openmpi/4.1.4 acme-netcdf/4.7.4/acme",
                  "export GATOR_INITIAL_MB=4000MB",
                  "export PATH=/ascldap/users/jgfouca/packages/valgrind-3.22.0/bin:$PATH",
                 ],
                 ["mpicxx","mpifort","mpicc"],
                  "",
                  "/sems-data-store/ACME/baselines/scream/master-baselines"),
    "lassen" : (["module --force purge", "module load git gcc/8.3.1 cuda/11.8.0 cmake/3.16.8 spectrum-mpi python/3.7.2", "export LLNL_USE_OMPI_VARS='y'",
                 "export PATH=/usr/gdata/e3sm/netcdf/bin:$PATH",
                 "export LD_LIBRARY_PATH=/usr/gdata/e3sm/netcdf/lib:$LD_LIBRARY_PATH",
                ],
                ["mpicxx","mpifort","mpicc"],
                 "bsub -Ip -qpdebug",
                 ""),
    "ruby-intel" : (["module --force purge", "module use --append /usr/gdata/e3sm/install/quartz/modulefiles", "module load StdEnv cmake/3.19.2 mkl/2022.1.0 intel-classic/2021.6.0-magic mvapich2/2.3.7 hdf5/1.12.2 netcdf-c/4.9.0 netcdf-fortran/4.6.0 parallel-netcdf/1.12.3 python/3.9.12 screamML-venv/0.0.1"],
                 ["mpicxx","mpifort","mpicc"],
                  "salloc --partition=pdebug",
                  ""),
    "quartz-intel" : (["module --force purge", "module use --append /usr/gdata/e3sm/install/quartz/modulefiles", "module load StdEnv cmake/3.19.2 mkl/2022.1.0 intel-classic/2021.6.0-magic mvapich2/2.3.7 hdf5/1.12.2 netcdf-c/4.9.0 netcdf-fortran/4.6.0 parallel-netcdf/1.12.3 python/3.9.12 screamML-venv/0.0.1"],
                 ["mpicxx","mpifort","mpicc"],
                  "salloc --partition=pdebug",
                  ""),
    "quartz-gcc" : (["module --force purge", "module load StdEnv cmake/3.16.8 mkl/2019.0 gcc/8.3.1 netcdf-fortran/4.4.4 netcdf/4.4.1.1 pnetcdf/1.9.0 mvapich2/2.3"],
                 ["mpicxx","mpifort","mpicc"],
                  "salloc --partition=pdebug",
                  ""),
    "syrah"  : (["module --force purge", "module load StdEnv cmake/3.16.8 mkl/2019.0 intel/19.0.4 netcdf-fortran/4.4.4 netcdf/4.4.1.1 pnetcdf/1.9.0 mvapich2/2.3"],
                 ["mpicxx","mpifort","mpicc"],
                  "salloc --partition=pbatch --time=60",
                  ""),
    "summit" : (["module purge", "module load cmake/3.18.4 gcc/9.3.0 spectrum-mpi/10.4.0.3-20210112 cuda/11.5.2 python/3.7-anaconda3 netcdf-c/4.8.0 netcdf-fortran/4.4.5 openblas/0.3.5","unset OMPI_CXX", "export OMPI_COMM_WORLD_RANK=0"],
                ["mpicxx","mpifort","mpicc"],
                "bsub -I -q batch -W 0:30 -P cli115 -nnodes 1",
                "/gpfs/alpine/cli115/proj-shared/scream/master-baselines"),
    "pm-cpu" : ([f"eval $({CIMEROOT}/CIME/Tools/get_case_env -c SMS.ne4pg2_ne4pg2.F2010-SCREAMv1.pm-cpu_gnu)"],
                ["CC","ftn","cc"],
                "salloc --time 00:30:00 --nodes=1 --constraint=cpu -q debug --account e3sm_g",
                "/global/cfs/cdirs/e3sm/baselines/gnu/scream/pm-cpu"),
    "pm-gpu" : ([f"eval $({CIMEROOT}/CIME/Tools/get_case_env -c SMS.ne4pg2_ne4pg2.F2010-SCREAMv1.pm-gpu_gnugpu)", "echo cuda=true"],
                ["CC","ftn","cc"],
                "salloc --time 02:00:00 --nodes=4 --constraint=gpu --gpus-per-node=4 --gpu-bind=none --exclusive -q regular --account e3sm_g",
                "/global/cfs/cdirs/e3sm/baselines/gnugpu/scream/pm-gpu"),
    "compy"   : (["module purge", "module load cmake/3.19.6 gcc/8.1.0  mvapich2/2.3.1 python/3.7.3"],
                 ["mpicxx","mpifort","mpicc"],
                  "srun --time 02:00:00 --nodes=1 -p short --exclusive --account e3sm",
                  ""),
    "chrysalis" : ([f"eval $({CIMEROOT}/CIME/Tools/get_case_env)", "export OMP_NUM_THREADS=1"],
                  ["mpic++","mpif90","mpicc"],
                  "srun --mpi=pmi2 -l -N 1 --kill-on-bad-exit --cpu_bind=cores",
                  "/lcrc/group/e3sm/baselines/chrys/intel/scream"),
    "anlgce"   : ([". /nfs/gce/software/spack/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/lmod-8.3-6fjdtku/lmod/lmod/init/sh", "module purge", "module load autoconf/2.69-bmnwajj automake/1.16.3-r7w24o4 libtool/2.4.6-uh3mpsu m4/1.4.19-7fztfyz cmake/3.20.5-zyz2eld gcc/11.1.0-qsjmpcg zlib/1.2.11-p7dmb5p", "export LD_LIBRARY_PATH=/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/mpich/4.0/gcc-11.1.0/lib:$LD_LIBRARY_PATH", "export PATH=/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/mpich/4.0/gcc-11.1.0/bin:/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/netcdf/4.8.0c-4.3.1cxx-4.5.3f-serial/gcc-11.1.0/bin:$PATH", "export NetCDF_ROOT=/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/netcdf/4.8.0c-4.3.1cxx-4.5.3f-serial/gcc-11.1.0", "export PERL5LIB=/nfs/gce/projects/climate/software/perl5/lib/perl5"],
                 ["mpicxx","mpifort","mpicc"],
                  "",
                  ""),
    "linux-generic" :        ([],["mpicxx","mpifort","mpicc"],"", ""),
    "linux-generic-debug" :  ([],["mpicxx","mpifort","mpicc"],"", ""),
    "linux-generic-serial" : ([],["mpicxx","mpifort","mpicc"],"", ""),
}

if pathlib.Path("~/.cime/scream_mach_specs.py").expanduser().is_file(): # pylint: disable=no-member
    sys.path.append(str(pathlib.Path("~/.cime").expanduser()))
    from scream_mach_specs import MACHINE_METADATA as LOCAL_MD # pylint: disable=import-error
    if len(LOCAL_MD) == 4:
        MACHINE_METADATA["local"] = LOCAL_MD
    else:
        print("WARNING! File '~/.cime/scream_mach_specs.py' was found, but the MACHINE_METADATA in there is badly formatted. Ignoring it.")

###############################################################################
def get_all_supported_machines():
###############################################################################
    return MACHINE_METADATA.keys()

###############################################################################
def is_machine_supported(machine):
###############################################################################
    return machine in MACHINE_METADATA.keys()

###############################################################################
def assert_machine_supported(machine):
###############################################################################
    expect(is_machine_supported(machine),
           "Machine {} is not currently supported by scream testing system.\n"
           " Note: you can also create a file `~/.cime/scream_mach_specs.py` with your local machine specs.".format(machine))

###############################################################################
def get_mach_env_setup_command(machine, ctest_j=None):
###############################################################################
    """
    ctest_j=None -> probe for hardware for testing resources
    ctest_j=-1   -> Skip CTEST_PARALLEL_LEVEL
    """
    assert_machine_supported(machine)

    mach_custom_env = MACHINE_METADATA[machine][0]
    if ctest_j != -1:
        ctest_j = get_mach_testing_resources(machine) if ctest_j is None else ctest_j
        mach_custom_env.append("export CTEST_PARALLEL_LEVEL={}".format(ctest_j))

    if not is_cuda_machine(machine):
        mach_custom_env.append("export OMP_PROC_BIND=spread")

    return mach_custom_env

###############################################################################
def get_mach_cxx_compiler(machine):
###############################################################################
    assert_machine_supported(machine)
    return MACHINE_METADATA[machine][1][0]

###############################################################################
def get_mach_f90_compiler(machine):
###############################################################################
    assert_machine_supported(machine)
    return MACHINE_METADATA[machine][1][1]

###############################################################################
def get_mach_c_compiler(machine):
###############################################################################
    assert_machine_supported(machine)
    return MACHINE_METADATA[machine][1][2]

###############################################################################
def get_mach_batch_command(machine):
###############################################################################
    assert_machine_supported(machine)
    return MACHINE_METADATA[machine][2]

###############################################################################
def get_mach_baseline_root_dir(machine):
###############################################################################
    """
    The pre-existing baselines root dir is used for integration testing only.
    It can be an empty string (""). If so, test-all-scream will use a default
    (ctest-build/baselines), and build baselines on the fly.
    This directory is interpreted as the directory where different builds
    subdirs are located (full_debug, full_sp_debug, debug_no_fpe).
    """
    assert_machine_supported(machine)
    return MACHINE_METADATA[machine][3]

###############################################################################
def logical_cores_per_physical_core():
###############################################################################
    return psutil.cpu_count() // psutil.cpu_count(logical=False)

###############################################################################
def get_available_cpu_count(logical=True):
###############################################################################
    """
    Get number of CPUs available to this process and its children. logical=True
    will include hyperthreads, logical=False will return only physical cores
    """
    affinity_len = len(psutil.Process().cpu_affinity())
    if not logical:
        hyperthread_ratio = logical_cores_per_physical_core()
        return int(affinity_len / hyperthread_ratio)
    else:
        return affinity_len

###############################################################################
def get_mach_compilation_resources():
###############################################################################
    return get_available_cpu_count()

###############################################################################
def get_mach_testing_resources(machine):
###############################################################################
    """
    The number of host cores is used to parallelize compilation,
    while the number of devices is used to parallelize testing.
    On CPU machines, the two will usually coincide, while on GPU
    machines they are going to be different (compile on CPU, run on GPU).
    One difference is that, for CPU machines, we allow hyperthreading for
    compilation but not for testing because we want to minimize fragmentation
    of jobs across cores.
    """
    if is_cuda_machine(machine):
        prefix = "srun " if is_salloc(machine) else ""
        return int(run_cmd_no_fail(f"{prefix}nvidia-smi -L | wc -l"))
    else:
        return get_available_cpu_count()

###############################################################################
def is_salloc(machine):
###############################################################################
    """
    Return true if we are running on an salloc'd job.
    """
    bcmd = get_mach_batch_command(machine)
    return "salloc" in bcmd and "srun" not in bcmd

###############################################################################
def is_cuda_machine(machine):
###############################################################################
    assert_machine_supported(machine)

    env_setup_raw = MACHINE_METADATA[machine][0]
    env_setup_str = " ".join(env_setup_raw)

    return "cuda" in env_setup_str.lower()

###############################################################################
def setup_mach_env(machine, ctest_j=None):
###############################################################################
    assert_machine_supported(machine)

    env_setup = get_mach_env_setup_command(machine, ctest_j=ctest_j)

    # Do something only if this machine has env specs
    if env_setup != []:
        # Running the env command only modifies the env in the subprocess
        # But we can return the resulting PATH, and update current env with that

        # Get the whole env string after running the env_setup command
        curr_env = run_cmd_no_fail("{{ {};  }} > /dev/null && env | sort".format(" && ".join(env_setup)),verbose=True)

        # Split by line. We are assuming that each env variable is *exactly* on one line
        curr_env_list = curr_env.split("\n")

        # For each line, split the string at the 1st '='.
        # The resulting length-2 stirng is (ENV_VAR_NAME, ENV_VAR_VALUE);
        # use it to update the os environment
        for item in curr_env_list:
            # On fedora systems, the environment contains the annoying entry (on 2 lines)
            #
            # BASH_FUNC_module()=() {  eval `/usr/bin/modulecmd bash $*`
            # }
            # Which breaks the assumption that each env var is on one line.
            # On some systems, this variable seems to have a different name,
            # and there can potentially be other BASH_FUNC_blah variables.
            # To get around this, discard lines that either do not contain '=',
            # or that start with BASH_FUNC_.
            if item.find("BASH_FUNC_") != -1 or item.find("=") == -1:
                continue

            # 2 means only 1st occurence will cause a split.
            # Just in case some env var value contains '='
            item_list = item.split("=",2)
            os.environ.update( dict( { item_list[0] : item_list[1] } ) )
