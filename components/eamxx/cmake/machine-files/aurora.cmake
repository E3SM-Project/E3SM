include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)
set(EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-np" CACHE STRING "" FORCE)
set(EKAT_MPI_EXTRA_ARGS "--label --cpu-bind depth -envall" CACHE STRING "")
set(EKAT_MPI_THREAD_FLAG "-d" CACHE STRING "")

set(CMAKE_Fortran_FLAGS " -fpscomp logicals -traceback -convert big_endian -assume byterecl -assume realloc_lhs -fp-model precise" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O2" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g" CACHE STRING "" FORCE)

SET(SYCL_COMPILE_FLAGS "-std=c++17 -fsycl -fsycl-device-code-split=per_kernel -fno-sycl-id-queries-fit-in-int -fsycl-unnamed-lambda")
SET(SYCL_LINK_FLAGS "-fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device pvc\" -flink-huge-device-code")

if (COMPILER MATCHES ".*gpu.*") # oneapi-ifxgpu
  include (${EKAT_MACH_FILES_PATH}/kokkos/intel-pvc.cmake)
  set(CMAKE_CXX_FLAGS  " -mlong-double-64 ${SYCL_COMPILE_FLAGS}" CACHE STRING "" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS  " -lifcore -mlong-double-64 ${SYCL_LINK_FLAGS} -fortlib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -fsycl-max-parallel-link-jobs=16 -Wl,--no-relax" CACHE STRING "" FORCE)
  set(SCREAM_MPI_ON_DEVICE ON CACHE STRING "" FORCE)
else() # oneapi-ifx
  include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  set(CMAKE_CXX_FLAGS  " -mlong-double-64" CACHE STRING "" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS  " -lifcore -mlong-double-64 -fortlib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" CACHE STRING "" FORCE)
endif()
