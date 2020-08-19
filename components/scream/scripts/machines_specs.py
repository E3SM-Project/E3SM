# MACHINE -> (env_setup, compiler, batch submit prefix, num host cores, num devices, pre-existing baselines root dir)

# Note: the number of host cores is used to parallelize compilation,
#       while the number of devices is used to parallelize testing.
#       On CPU machines, the two will usually coincide, while on GPU
#       machines they are going to be different (compile on CPU, run on GPU).

# Note: the pre-existing baselines root dir is used for integration testing only.
#       It can be an empty string (""). If so, test-all-scream will use a default
#       (ctest-build/baselines), and build baselines on the fly.
#       This directory is interpreted as the directory where different builds
#       subdirs are located (full_debug, full_sp_debug, debug_no_fpe).

from utils import expect, get_cpu_core_count, run_cmd_no_fail
import os

MACHINE_METADATA = {
    "melvin"   : (["module purge", "module load sems-env", "module load sems-gcc/7.3.0 sems-openmpi/1.10.1 sems-gcc/7.3.0 sems-git/2.10.1 sems-cmake/3.12.2 sems-python/3.5.2"],
                  "$(which mpicxx)",
                  "",
                  24,
                  24,
                  ""),
    "bowman"   : (["module purge", "module load openmpi/1.10.6/intel/17.2.174 git/2.8.2 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-bowman/bin:$PATH"],
                  "$(which mpicxx)",
                  "srun",
                  68,
                  68,
                  "/home/projects/e3sm/scream/pr-autotester/master-baselines/bowman/"),
    "blake"    : (["module purge", "module load openmpi/2.1.5/intel/19.1.144 git/2.9.4 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-blake/bin:$PATH"],
                  "$(which mpicxx)",
                  "srun",
                  48,
                  48,
                  "/home/projects/e3sm/scream/pr-autotester/master-baselines/blake/"),
    "weaver"   : (["module purge", "module load devpack/20190814/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105 git/2.10.1 python/3.7.3", "module switch cmake/3.18.0"],
                  "$(which mpicxx)",
                  "bsub -I -q rhel7W",
                  40,
                  4,
                  "/home/projects/e3sm/scream/pr-autotester/master-baselines/weaver/"),
    "white"    : (["module purge", "module load devpack/20181011/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88 git/2.10.1", "module switch cmake/3.12.3 cmake/3.18.0", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-white/bin:$PATH"],
                  "$(which mpicxx)",
                  "bsub -I -q rhel7G",
                  64,
                  4,
                  "/home/projects/e3sm/scream/pr-autotester/master-baselines/white/"),
    "lassen" : (["module --force purge", "module load git gcc/7.3.1 cuda/10.1.243 cmake/3.14.5 spectrum-mpi netcdf/4.7.0 python/3.7.2", "export LLNL_USE_OMPI_VARS='y'"],
                  "$(which mpicxx)",
                  "bsub -Ip",
                  44,
                  4,
                  ""),
    "quartz" : (["module --force purge", "module load StdEnv cmake/3.14.5 mkl/2019.0 intel/19.0.4 netcdf-fortran/4.4.4 netcdf/4.4.1.1"],
                  "$(which mpicxx)",
                  "salloc --partition=pdebug",
                  36,
                  36,
                  ""),
    "syrah"  : (["module --force purge", "module load StdEnv cmake/3.14.5 mkl/2019.0 intel/19.0.4 netcdf-fortran/4.4.4 netcdf/4.4.1.1"],
                  "$(which mpicxx)",
                  "salloc --partition=pdebug",
                  16,
                  16,
                  ""),
    "summit" : (["module purge", "module load cmake/3.15.2 gcc/6.4.0 spectrum-mpi/10.3.0.1-20190611 cuda/10.1.168 python/3.6.6-anaconda3-5.3.0"],
                "$(which mpicxx)",
                "bsub -I -q batch -W 0:30 -P cli115 -nnodes 1",
                44,
                6,
                ""),
    "cori"   : (["eval $(../../cime/scripts/Tools/get_case_env)", "export OMP_NUM_THREADS=68"],
                "$(which CC)",
                "srun --time 02:00:00 --nodes=1 --constraint=knl,quad,cache --exclusive -q regular --account e3sm",
                68,
                68,
                ""),
    "compy"   : (["module purge", "module load cmake/3.11.4 gcc/8.1.0  mvapich2/2.3.1 python/3.7.3"],
                  "$(which mpicxx)",
                  "",
                  40,
                  40,
                  ""),

    "linux-generic" : ([],"$(which mpicxx)","", get_cpu_core_count(), get_cpu_core_count(),""),
    "linux-generic-debug" : ([],"$(which mpicxx)","", get_cpu_core_count(), get_cpu_core_count(),""),
    "linux-generic-serial" : ([],"$(which mpicxx)","", get_cpu_core_count(), get_cpu_core_count(),""),
}

###############################################################################
def is_machine_supported (machine):
###############################################################################

    return machine in MACHINE_METADATA.keys()

###############################################################################
def get_mach_env_setup_command (machine):
###############################################################################

    expect (is_machine_supported(machine), "Error! Machine {} is not currently supported by scream testing system.".format(machine))

    return MACHINE_METADATA[machine][0]

###############################################################################
def get_mach_cxx_compiler (machine):
###############################################################################

    expect (is_machine_supported(machine), "Error! Machine {} is not currently supported by scream testing system.".format(machine))

    return MACHINE_METADATA[machine][1]

###############################################################################
def get_mach_batch_command (machine):
###############################################################################

    expect (is_machine_supported(machine), "Error! Machine {} is not currently supported by scream testing system.".format(machine))

    return MACHINE_METADATA[machine][2]

###############################################################################
def get_mach_compilation_resources (machine):
###############################################################################

    expect (is_machine_supported(machine), "Error! Machine {} is not currently supported by scream testing system.".format(machine))

    return MACHINE_METADATA[machine][3]

###############################################################################
def get_mach_testing_resources (machine):
###############################################################################

    expect (is_machine_supported(machine), "Error! Machine {} is not currently supported by scream testing system.".format(machine))

    return MACHINE_METADATA[machine][4]


###############################################################################
def get_mach_baseline_root_dir (machine,default_dir):
###############################################################################

    expect (is_machine_supported(machine), "Error! Machine {} is not currently supported by scream testing system.".format(machine))

    if MACHINE_METADATA[machine][5]=="":
        return default_dir
    else:
        return MACHINE_METADATA[machine][5]

###############################################################################
def setup_mach_env (machine):
###############################################################################

    expect (is_machine_supported(machine), "Error! Machine {} is not currently supported by scream testing system.".format(machine))

    env_setup = get_mach_env_setup_command(machine)

    # Do something only if this machine has env specs
    if env_setup != []:
        # Running the env command only modifies the env in the subprocess
        # But we can return the resulting PATH, and update current env with that

        # Get the whole env string after running the env_setup command
        curr_env = run_cmd_no_fail("{{ {};  }} > /dev/null && env | sort".format(";".join(env_setup)))

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
            # To get around this, discard lines that either do not contain '=',
            # or that start with BASH_FUNC_module.
            if item.find("BASH_FUNC_module") != -1 or item.find("=") == -1:
                continue


            # 2 means only 1st occurence will cause a split.
            # Just in case some env var value contains '='
            item_list = item.split("=",2)
            os.environ.update( dict( { item_list[0] : item_list[1] } ) )
