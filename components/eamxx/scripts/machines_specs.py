from utils import expect, run_cmd_no_fail, ensure_psutil

import os, sys, pathlib
ensure_psutil()
import psutil

EAMXX_DIR = pathlib.Path(__file__).parent.parent
CIMEROOT = os.path.join(EAMXX_DIR,"..","..","cime")

###############################################################################
def logical_cores_per_physical_core():
###############################################################################
    return psutil.cpu_count() // psutil.cpu_count(logical=False)

###############################################################################
def get_cpu_ids_from_slurm_env_var():
###############################################################################
    """
    Parse the SLURM_CPU_BIND_LIST, and use the hexadecimal value to determine
    which CPUs on this node are assigned to the job
    NOTE: user should check that the var is set BEFORE calling this function
    """

    cpu_bind_list = os.getenv('SLURM_CPU_BIND_LIST')

    expect (cpu_bind_list is not None,
            "SLURM_CPU_BIND_LIST environment variable is not set. Check, before calling this function")

    # Remove the '0x' prefix and convert to an integer
    mask_int = int(cpu_bind_list, 16)

    # Generate the list of CPU IDs
    cpu_ids = []
    for i in range(mask_int.bit_length()):  # Check each bit position
        if mask_int & (1 << i):  # Check if the i-th bit is set
            cpu_ids.append(i)

    return cpu_ids

###############################################################################
def get_available_cpu_count(logical=True):
###############################################################################
    """
    Get number of CPUs available to this process and its children. logical=True
    will include hyperthreads, logical=False will return only physical cores
    """
    if 'SLURM_CPU_BIND_LIST' in os.environ:
        cpu_count = len(get_cpu_ids_from_slurm_env_var())
    else:
        cpu_count = len(psutil.Process().cpu_affinity())

    if not logical:
        hyperthread_ratio = logical_cores_per_physical_core()
        return int(cpu_count / hyperthread_ratio)
    else:
        return cpu_count

###############################################################################
class Machine(object):
###############################################################################
    """
    Parent class for objects describing a machine to use for EAMxx standalone testing.
    """
    concrete = False
    name = ""
    mach_file = ""
    env_setup = []
    gpu_arch = "none"
    num_bld_res = -1
    num_run_res = -1
    batch = ""
    cxx_compiler = "mpicxx"
    c_compiler   = "mpicc"
    ftn_compiler = "mpifort"
    baselines_dir = ""

    @classmethod
    def setup_base(cls,name,num_bld_res=-1,num_run_res=-1):
        cls.name = name
        cls.mach_file = pathlib.Path(EAMXX_DIR) / "cmake" / "machine-files" / (cls.name + ".cmake")
        cls.num_bld_res = num_bld_res if num_bld_res > 0 else get_available_cpu_count()
        cls.num_run_res = num_run_res if num_run_res > 0 else get_available_cpu_count()

    @classmethod
    def uses_gpu (cls):
        return cls.gpu_arch!="none"

###############################################################################
class Generic(Machine):
###############################################################################
    concrete = True

    @classmethod
    def setup(cls):
        super().setup_base("linux-generic")
        # Set baselines_dir to PWD/../ctest-build/baselines
        cls.baselines_dir = os.path.join(os.path.dirname(__file__), '..', 'ctest-build', 'baselines')

###############################################################################
class CrayMachine(Machine):
###############################################################################
    @classmethod
    def setup_cray(cls,name):
        super().setup_base(name)

        cls.cxx_compiler = "CC"
        cls.c_compiler   = "cc"
        cls.ftn_compiler = "ftn"

###############################################################################
class PM(CrayMachine):
###############################################################################
    @classmethod
    def setup_pm(cls,partition):
        expect (partition in ['cpu', 'gpu'], "Unknown Perlmutter partition")

        super().setup_cray("pm-"+partition)

        compiler = "gnu" if partition=="cpu" else "gnugpu"

        cls.env_setup = [f"eval $({CIMEROOT}/CIME/Tools/get_case_env -c SMS.ne4pg2_ne4pg2.F2010-SCREAMv1.{cls.name}_{compiler})"]
        cls.batch = f"salloc --account e3sm_g --constraint={partition}"
        if partition=="cpu":
            cls.batch += "--time 00:30:00 --nodes=1 -q debug"
        else:
            cls.batch += "--time 02:00:00 --nodes=4 --gpus-per-node=4 --gpu-bind=none --exclusive -q regular"
        
        cls.baselines_dir = f"/global/cfs/cdirs/e3sm/baselines/{compiler}/scream/{cls.name}"

###############################################################################
class PMCPU(PM):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_pm("cpu")

###############################################################################
class PMGPU(PM):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_pm("gpu")

        cls.num_run_res = 4 # four gpus
        cls.gpu_arch = "cuda"

###############################################################################
class Chrysalis(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("chrysalis")

        cls.env_setup = [f"eval $({CIMEROOT}/CIME/Tools/get_case_env)", "export OMP_NUM_THREADS=1"]
        cls.batch = "srun --mpi=pmi2 -l -N 1 --kill-on-bad-exit --cpu_bind=cores"
        cls.baselines_dir = "/lcrc/group/e3sm/baselines/chrys/intel/scream"

        cls.cxx_compiler = "mpic++"
        cls.c_compiler   = "mpicc"
        cls.ftn_compiler = "mpif90"

###############################################################################
class Mappy(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("mappy")

        cls.env_setup = ["module purge",
                          "module load sems-cmake/3.27.9 sems-git/2.42.0 sems-gcc/11.4.0 sems-openmpi-no-cuda/4.1.6 sems-netcdf-c/4.9.2 sems-netcdf-cxx/4.2 sems-netcdf-fortran/4.6.1 sems-parallel-netcdf/1.12.3 sems-openblas",
                          "export GATOR_INITIAL_MB=4000MB",
                         ]
        cls.baselines_dir = "/sems-data-store/ACME/baselines/scream/master-baselines"

###############################################################################
class Weaver(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("weaver")

        cls.env_setup = ["source /etc/profile.d/modules.sh",
                          "module purge",
                          "module load cmake/3.25.1 git/2.39.1 python/3.10.8 py-netcdf4/1.5.8 gcc/11.3.0 cuda/11.8.0 openmpi netcdf-c netcdf-fortran parallel-netcdf netlib-lapack",
                          "export HDF5_USE_FILE_LOCKING=FALSE"
                         ]
        cls.baselines_dir = "/home/projects/e3sm/scream/pr-autotester/master-baselines/weaver/"
        cls.batch = "bsub -I -q rhel8 -n 4 -gpu num=4"

        cls.num_run_res = 4 # four gpus
        cls.gpu_arch = "cuda"

###############################################################################
class Compy(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("compy")

        cls.env_setup = ["module purge", "module load cmake/3.19.6 gcc/8.1.0  mvapich2/2.3.1 python/3.7.3"]
        cls.batch = "srun --time 02:00:00 --nodes=1 -p short --exclusive --account e3sm"

###############################################################################
class GHCISNLCPU(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("ghci-snl-cpu")
        cls.baselines_dir = "/projects/e3sm/baselines/scream/ghci-snl-cpu"

###############################################################################
class GHCISNLCuda(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base(name="ghci-snl-cuda",num_bld_res=16,num_run_res=1)
        cls.baselines_dir = "/projects/e3sm/baselines/scream/ghci-snl-cuda"
        cls.gpu_arch = "cuda"

###############################################################################
class GHCIOCI(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base(name="ghci-oci")
        cls.env_setup = [f"eval $({CIMEROOT}/CIME/Tools/get_case_env -c SMS.ne4pg2_ne4pg2.F2010-SCREAMv1.ghci-oci_gnu)"]

###############################################################################
class Lassen(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("lassen")
        cls.baselines_dir = "/projects/e3sm/baselines/scream/master-baselines"

        cls.env_setup = ["module --force purge",
                          "module load git gcc/8.3.1 cuda/11.8.0 cmake/3.16.8 spectrum-mpi python/3.7.2",
                          "export LLNL_USE_OMPI_VARS='y'",
                          "export PATH=/usr/workspace/e3sm/netcdf/bin:$PATH",
                          "export LD_LIBRARY_PATH=/usr/workspace/e3sm/netcdf/lib:$LD_LIBRARY_PATH",
                         ]
        cls.batch = "bsub -Ip -qpdebug"

        cls.num_run_res = 4 # four gpus
        cls.gpu_arch = "cuda"

###############################################################################
class LLNLIntel(Machine):
###############################################################################
    @classmethod
    def setup_llnl_intel(cls,name):
        super().setup_base(name)

        cls.env_setup = ["module --force purge",
                          "module use --append /usr/workspace/e3sm/install/quartz/modulefiles",
                          "module load StdEnv cmake/3.19.2 mkl/2022.1.0 intel-classic/2021.6.0-magic mvapich2/2.3.7 hdf5/1.12.2 netcdf-c/4.9.0 netcdf-fortran/4.6.0 parallel-netcdf/1.12.3 python/3.9.12 screamML-venv/0.0.1"
                         ]
        cls.batch = "salloc --partition=pdebug",


###############################################################################
class RubyIntel(LLNLIntel):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_llnl_intel("ruby-intel")

###############################################################################
class DaneIntel(LLNLIntel):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_llnl_intel("dane-intel")

###############################################################################
class QuartzIntel(LLNLIntel):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_llnl_intel("quartz-intel")

        cls.env_setup = ["module --force purge",
                         "module use --append /usr/workspace/e3sm/install/quartz/modulefiles",
                         "module load StdEnv cmake/3.19.2 mkl/2022.1.0 intel-classic/2021.6.0-magic mvapich2/2.3.7 hdf5/1.12.2 netcdf-c/4.9.0 netcdf-fortran/4.6.0 parallel-netcdf/1.12.3 python/3.9.12 screamML-venv/0.0.1"
                        ]
        cls.batch = "salloc --partition=pdebug",

###############################################################################
class QuartzGCC(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("quartz-gcc")

        cls.env_setup = ["module --force purge",
                         "module load StdEnv cmake/3.16.8 mkl/2019.0 gcc-8.3.1 netcdf-fortran/4.4.4 netcdf/4.4.1.1 pnetcdf/1.9.0 mvapich2/2.3"
                        ]
        cls.batch = "salloc --partition=pdebug",

###############################################################################
class Syrah(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("syrah")

        cls.env_setup = ["module --force purge",
                         "module load StdEnv cmake/3.16.8 mkl/2019.0 intel/19.0.4 netcdf-fortran/4.4.4 netcdf/4.4.1.1 pnetcdf/1.9.0 mvapich2/2.3"
                        ]
        cls.batch = "salloc --partition=pdebug --time=60",

###############################################################################
class AnlGce(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("anlgce")

        cls.env_setup = [". /nfs/gce/software/spack/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/lmod-8.3-6fjdtku/lmod/lmod/init/sh",
                         "module purge",
                         "module load autoconf/2.69-bmnwajj automake/1.16.3-r7w24o4 libtool/2.4.6-uh3mpsu m4/1.4.19-7fztfyz cmake/3.20.5-zyz2eld gcc/11.1.0-qsjmpcg zlib/1.2.11-p7dmb5p",
                         "export LD_LIBRARY_PATH=/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/mpich/4.0/gcc-11.1.0/lib:$LD_LIBRARY_PATH",
                         "export PATH=/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/mpich/4.0/gcc-11.1.0/bin:/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/netcdf/4.8.0c-4.3.1cxx-4.5.3f-serial/gcc-11.1.0/bin:$PATH",
                         "export NetCDF_ROOT=/nfs/gce/projects/climate/software/linux-ubuntu20.04-x86_64/netcdf/4.8.0c-4.3.1cxx-4.5.3f-serial/gcc-11.1.0",
                         "export PERL5LIB=/nfs/gce/projects/climate/software/perl5/lib/perl5"
                        ]

###############################################################################
class AnlGceUb22(Machine):
###############################################################################
    concrete = True
    @classmethod
    def setup(cls):
        super().setup_base("anlgce-ub22")

        cls.env_setup = [". /nfs/gce/software/custom/linux-ubuntu22.04-x86_64/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-11.2.0/lmod-8.5.6-hkjjxhp/lmod/lmod/init/sh",
                         "module purge",
                         "module load gcc/12.1.0",
                         "export LD_LIBRARY_PATH=/nfs/gce/projects/climate/software/linux-ubuntu22.04-x86_64/mpich/4.1.2/gcc-12.1.0/lib:$LD_LIBRARY_PATH",
                         "export PATH=/nfs/gce/projects/climate/software/linux-ubuntu22.04-x86_64/mpich/4.1.2/gcc-12.1.0/bin:/nfs/gce/projects/climate/software/linux-ubuntu22.04-x86_64/netcdf/4.8.0c-4.3.1cxx-4.5.3f-serial/gcc-12.1.0/bin:$PATH",
                         "export NetCDF_ROOT=/nfs/gce/projects/climate/software/linux-ubuntu22.04-x86_64/netcdf/4.8.0c-4.3.1cxx-4.5.3f-serial/gcc-12.1.0",
                         "export PERL5LIB=/nfs/gce/projects/climate/software/perl5/lib/perl5"
                        ]

###############################################################################
def get_all_machines (base=Machine):
###############################################################################
    # If the user has the file ~/.cime/scream_mach_specs.py, import the machine type Local from it
    if pathlib.Path("~/.cime/scream_mach_specs.py").expanduser().is_file(): # pylint: disable=no-member
        sys.path.append(str(pathlib.Path("~/.cime").expanduser()))
        # Add the directory of this file to sys.path, so that ~/.cime/scream_mach_specs.py,
        # if present, can simply do "from machine_specs import Machine" to define machine Local
        sys.path.append(os.path.dirname(__file__))

        # Pylint does not see that this import adds one subclass to Machine
        # Also, its static analysis runs on a static sys.path, without ~/.cime,
        # so it doesn't find the module
        from scream_mach_specs import Local# pylint: disable=unused-import, import-error

    machines = []
    for m in base.__subclasses__():
        if m.concrete:
            m.setup() # Init the class static data
            machines.append(m)

        machines.extend(get_all_machines(m))

    return machines

###############################################################################
def get_machine (name):
###############################################################################

    all_machines = get_all_machines()

    for m in all_machines:
        if m.name==name:
            return m

    raise ValueError(f"Machine with name '{name}' not found.")

###############################################################################
def get_all_supported_machines ():
###############################################################################
    return [m.name for m in get_all_machines()]

###############################################################################
def is_machine_supported(machine):
###############################################################################
    return machine in get_all_supported_machines()

###############################################################################
def assert_machine_supported(machine):
###############################################################################
    expect(is_machine_supported(machine),
           "Machine {} is not currently supported by scream testing system.\n"
          f" Currently supported machines: {','.join(get_all_supported_machines())}\n"
           " Note: you can also create a file `~/.cime/scream_mach_specs.py` with your local machine specs.".format(machine))

###############################################################################
def get_mach_env_setup_command(machine, ctest_j=None):
###############################################################################
    """
    ctest_j=None -> probe for hardware for testing resources
    ctest_j=-1   -> Skip CTEST_PARALLEL_LEVEL
    """

    mach = get_machine(machine)
    mach_custom_env = mach.env_setup
    if ctest_j != -1:
        ctest_j = mach.num_run_res if ctest_j is None else ctest_j
        mach_custom_env.append("export CTEST_PARALLEL_LEVEL={}".format(ctest_j))

    if not mach.uses_gpu():
        mach_custom_env.append("export OMP_PROC_BIND=spread")

    return mach_custom_env

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
