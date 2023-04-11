string(APPEND FFLAGS " -Wno-implicit-interface")

if (NOT DEBUG)
  string(APPEND FFLAGS " -O2")
  string(APPEND CFLAGS " -O2")
endif()
string(APPEND CXX_LIBS " -lstdc++")

string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")
if (NOT MPILIB STREQUAL mpi-serial)
  string(APPEND SLIBS " -L$ENV{ADIOS2_DIR}/lib64 -ladios2_c_mpi -ladios2_c -ladios2_core_mpi -ladios2_core -ladios2_evpath -ladios2_ffs -ladios2_dill -ladios2_atl -ladios2_enet")
endif()
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
