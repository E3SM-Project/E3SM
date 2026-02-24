include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)
set(EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-np" CACHE STRING "" FORCE)
set(EKAT_MPI_EXTRA_ARGS "--label --cpu-bind depth -envall" CACHE STRING "")
set(EKAT_MPI_THREAD_FLAG "-d" CACHE STRING "")

SET(SYCL_COMPILE_FLAGS "-std=c++17 -fsycl -fsycl-device-code-split=per_kernel -fno-sycl-id-queries-fit-in-int -fsycl-unnamed-lambda")
SET(SYCL_LINK_FLAGS "-fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\"")

if (COMPILER MATCHES ".*gpu.*") # oneapi-ifxgpu
  include (${EKAT_MACH_FILES_PATH}/kokkos/intel-pvc.cmake)
  set(CMAKE_CXX_FLAGS  " --intel -mlong-double-64 ${SYCL_COMPILE_FLAGS}" CACHE STRING "" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS  " -lifcore --intel -lsycl -mlong-double-64 ${SYCL_LINK_FLAGS} -fortlib" CACHE STRING "" FORCE)
else() # oneapi-ifx
  include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  set(CMAKE_CXX_FLAGS  " --intel -mlong-double-64" CACHE STRING "" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS  " -lifcore --intel -mlong-double-64 -fortlib" CACHE STRING "" FORCE)
endif()
