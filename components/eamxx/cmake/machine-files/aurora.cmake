include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/intel-pvc.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)
set(EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-np" CACHE STRING "" FORCE)
set(EKAT_MPI_EXTRA_ARGS "--label --cpu-bind depth -envall" CACHE STRING "")
set(EKAT_MPI_THREAD_FLAG "-d" CACHE STRING "")

SET(SYCL_COMPILE_FLAGS "-std=c++17 -fsycl -fsycl-device-code-split=per_kernel -fno-sycl-id-queries-fit-in-int -fsycl-unnamed-lambda")
SET(SYCL_LINK_FLAGS "-fsycl -fsycl-link-huge-device-code -fsycl-device-code-split=per_kernel  -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\"")

set(CMAKE_CXX_FLAGS " -\-intel -Xclang -fsycl-allow-virtual-functions -mlong-double-64 -O3 -DNDEBUG ${SYCL_COMPILE_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "-fc=ifx -O3 -DNDEBUG -DCPRINTEL -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS " -lifcore -\-intel -Xclang -fsycl-allow-virtual-functions -lsycl -mlong-double-64 -DNDEBUG ${SYCL_LINK_FLAGS} -fortlib" CACHE STRING "" FORCE)

#this is needed for cime builds!
set(NETCDF_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2024.04.15.002")
set(NETCDF_DIR  "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2024.04.15.002")
set(NETCDF_C_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2024.04.15.002")
set(NETCDF_C      "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2024.04.15.002")
#this one is for rrtmgp
set(NetCDF_C_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2024.04.15.002" CACHE STRING "")
set(NETCDF_FORTRAN_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2024.04.15.002")
set(PNETCDF_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/pnetcdf/1.12.3/oneapi.eng.2024.04.15.002")


