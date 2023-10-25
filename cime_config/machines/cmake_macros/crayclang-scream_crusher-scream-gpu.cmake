if (compile_threaded)
  string(APPEND CFLAGS " -fopenmp")
  string(APPEND FFLAGS " -fopenmp")
  string(APPEND CXXFLAGS " -fopenmp")
  string(APPEND LDFLAGS " -fopenmp")
endif()

set(PIO_FILESYSTEM_HINTS "gpfs")

set(MPICXX "hipcc")
set(SCXX "hipcc")

string(APPEND CXXFLAGS " -I${MPICH_DIR}/include")
string(APPEND LDFLAGS " -L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa")

# For YAKL's -lroctx64 -lrocfft; the rocm module doesn't set this.
string(APPEND LDFLAGS " -L$ENV{ROCM_PATH}/lib")

# 'NOT DEBUG': this resolves a crash in mct in docn init
# 'DEBUG' casee, too: resolves a build error in elm/src/main/elm_varctl.F90 due to several OpenACC syntax errors
#if (NOT DEBUG)
  string(APPEND CFLAGS " -O2 -hnoacc -hfp0 -hipa0")
  string(APPEND FFLAGS " -O2 -hnoacc -hfp0 -hipa0")
#endif()

string(APPEND CPPDEFS " -DCPRCRAY")

#set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "See SCREAM issue #2080.")
