from machines_specs import Machine

###############################################################################
class Local(Machine):
###############################################################################
  concrete = True

  name = "local"
  os_type = "macos"
  num_bld_res = 8
  num_run_res = 8

  @classmethod
  def setup(cls):
    super().setup_base()

    import socket
    my_host = socket.gethostname()
    print('==========================')
    print(f'my_host = {my_host}')
    print('==========================')

    if "s1088599" in my_host:
      cls.env_setup = ["source ${HOME}/.cime/cpu/scream-setup.sh"]
      cls.mach_file = "${HOME}/.cime/cpu/scream_mach_file.cmake"
      print('==========================')
      print(f'mach_file path = {cls.mach_file}')
      print('==========================')
    else:
      cls.env_setup = ["source ${HOME}/.cime/cpu/scream-setup.sh"]
      cls.mach_file = "${HOME}/.cime/cpu/scream_mach_file.cmake"

    # cls.baselines_dir = "/ascldap/users/mjschm/scream-data/master-baselines"


    # tpl_path="/home/mjschm/TPLs"
    # mpi_bindir = f"/projects/aue/cee/deploy/e1059ca4/linux-rhel8-x86_64/gcc-11.4.0/openmpi-4.1.6-dl5lbds/bin"
    # cls.cxx_compiler = "g++-15"
    # cls.cxx_compiler = "clang++"
    cls.cxx_compiler = "mpicxx"
    # cls.c_compiler = "gcc-15"
    # cls.c_compiler = "clang"
    cls.c_compiler = "mpicc"
    # cls.ftn_compiler = "gfortran"
    cls.ftn_compiler = "mpif90"

    # cls.gpu_arch = "cuda"


# MACHINE_METADATA = (
#     ["source ${HOME}/.cime/scream-setup.sh"],
#     ["mpicxx","mpifort","mpicc"],
#     "",
#     "/ascldap/users/mjschm/scream-data/master-baselines"
# )
