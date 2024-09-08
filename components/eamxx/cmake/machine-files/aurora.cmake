include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/intel-pvc.cmake)
# kokkos sycl is on in the above file
#include (${EKAT_MACH_FILES_PATH}/kokkos/sycl.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

#AB flags from ekat
# -fsycl -fsycl-unnamed-lambda -sycl-std=2020 -qopenmp-simd -Wsycl-strict -fsycl-device-code-split=per_kernel
SET(SYCL_COMPILE_FLAGS "-std=c++17 -fsycl -fsycl-device-code-split=per_kernel -fno-sycl-id-queries-fit-in-int -fsycl-unnamed-lambda")
SET(SYCL_LINK_FLAGS "-fsycl -fsycl-link-huge-device-code -fsycl-device-code-split=per_kernel  -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\"")

#SET(MPICH_DIR "/soft/restricted/CNDA/updates/mpich/52.2/mpich-ofi-all-icc-default-pmix-gpu-drop52/" CACHE STRING "")

set(CMAKE_CXX_FLAGS " -\-intel -Xclang -fsycl-allow-virtual-functions -mlong-double-64 -O3 -DNDEBUG ${SYCL_COMPILE_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "-fc=ifx -O3 -DNDEBUG -DCPRINTEL -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS " -lifcore -\-intel -Xclang -fsycl-allow-virtual-functions -lsycl -mlong-double-64 -DNDEBUG ${SYCL_LINK_FLAGS} -fortlib" CACHE STRING "" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS " -Wl,-\-defsym,main=MAIN_\_ -lifcore -\-intel -Xclang -fsycl-allow-virtual-functions -lsycl -mlong-double-64 -DNDEBUG ${SYCL_LINK_FLAGS} -fortlib  -L${MPICH_DIR}/lib" CACHE STRING "" FORCE)



#this is needed for cime builds!
set(NETCDF_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2024.04.15.002")
set(NETCDF_DIR  "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2024.04.15.002")
set(NETCDF_C_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2024.04.15.002")
set(NETCDF_C      "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2024.04.15.002")
#this one is for rrtmgp
set(NetCDF_C_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2024.04.15.002" CACHE STRING "")
set(NETCDF_FORTRAN_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2024.04.15.002")
set(PNETCDF_PATH "/lus/flare/projects/CSC249ADSE15_CNDA/software/pnetcdf/1.12.3/oneapi.eng.2024.04.15.002")


