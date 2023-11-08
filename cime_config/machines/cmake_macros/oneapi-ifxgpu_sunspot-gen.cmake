
set(CXX_LINKER "CXX")

execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nf-config --flibs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)

string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0} -Wl,-rpath -Wl,$ENV{NETCDF_PATH}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nc-config --libs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)

string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0}")
string(APPEND SLIBS " -fiopenmp -fopenmp-targets=spir64")

set(NETCDF_PATH "$ENV{NETCDF_PATH}")
set(PNETCDF_PATH "$ENV{PNETCDF_PATH}")

set(USE_SYCL "TRUE")

string(APPEND KOKKOS_OPTIONS " -DCMAKE_CXX_STANDARD=17 -DKokkos_ENABLE_SERIAL=On -DKokkos_ARCH_INTEL_GEN=On -DKokkos_ENABLE_SYCL=On -DKokkos_ENABLE_EXPLICIT_INSTANTIATION=Off")

string(APPEND SYCL_FLAGS " -\-intel -Xclang -fsycl-allow-virtual-functions -fsycl -mlong-double-64 -fsycl-device-code-split=per_kernel -fno-sycl-id-queries-fit-in-int -fsycl-unnamed-lambda")

#string(APPEND SYCL_FLAGS " -\-intel -fsycl")
string(APPEND CXX_LDFLAGS " -Wl,-\-defsym,main=MAIN_\_ -lifcore -\-intel -Xclang -fsycl-allow-virtual-functions -fsycl -lsycl -mlong-double-64 -fsycl-link-huge-device-code -fsycl-device-code-split=per_kernel -fsycl-targets=spir64")

SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "")
SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_FORTRAN_COMPILER "mpifort" CACHE STRING "")



