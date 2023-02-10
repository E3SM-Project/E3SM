if (COMP_NAME STREQUAL elm)
  # See Land NaNs in conditionals: https://github.com/E3SM-Project/E3SM/issues/4996
  string(APPEND FFLAGS " -hfp0")
endif()
# -em -ef generates modulename.mod (lowercase files) to support
# Scorpio installs
# Disable ipa and zero initialization are for other NaN isues:
# https://github.com/E3SM-Project/E3SM/pull/5208
string(APPEND FFLAGS " -hipa0 -hzero -em -ef -hnoacc")

set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
string(APPEND CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")
string(APPEND CXX_LIBS " -lstdc++")

string(APPEND CXXFLAGS " -I$ENV{MPICH_DIR}/include --amdgpu-target=gfx90a")
string(APPEND SLIBS    " -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

set(USE_HIP "TRUE")
string(APPEND HIP_FLAGS "-munsafe-fp-atomics -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=$ENV{ROCM_PATH} --offload-arch=gfx90a -x hip")
