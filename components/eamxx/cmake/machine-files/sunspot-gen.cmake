include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/intel-gen.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/sycl.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

#AB flags from ekat
# -fsycl -fsycl-unnamed-lambda -sycl-std=2020 -qopenmp-simd -Wsycl-strict -fsycl-device-code-split=per_kernel
SET(SYCL_COMPILE_FLAGS "-std=c++17 -fsycl -fsycl-device-code-split=per_kernel -fno-sycl-id-queries-fit-in-int -fsycl-unnamed-lambda")
SET(SYCL_LINK_FLAGS "-fsycl -fsycl-link-huge-device-code -fsycl-device-code-split=per_kernel -fsycl-targets=spir64")

#SET(MPICH_DIR "/soft/restricted/CNDA/updates/mpich/52.2/mpich-ofi-all-icc-default-pmix-gpu-drop52/" CACHE STRING "")

set(CMAKE_CXX_FLAGS " -\-intel -Xclang -fsycl-allow-virtual-functions -mlong-double-64 -O3 -DNDEBUG ${SYCL_COMPILE_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "-fc=ifx -O3 -DNDEBUG -DCPRINTEL -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS " -lifcore -\-intel -Xclang -fsycl-allow-virtual-functions -lsycl -mlong-double-64 -DNDEBUG ${SYCL_LINK_FLAGS} -fortlib" CACHE STRING "" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS " -Wl,-\-defsym,main=MAIN_\_ -lifcore -\-intel -Xclang -fsycl-allow-virtual-functions -lsycl -mlong-double-64 -DNDEBUG ${SYCL_LINK_FLAGS} -fortlib  -L${MPICH_DIR}/lib" CACHE STRING "" FORCE)



set(NETCDF_PATH "$ENV{NETCDF_PATH}")
set(NETCDF_C_PATH "$ENV{NETCDF_PATH}")
#this one is for rrtmgp
set(NetCDF_C_PATH "$ENV{NETCDF_PATH}" CACHE STRING "")
set(NETCDF_FORTRAN_PATH "$ENV{NETCDF_PATH}")
set(PNETCDF_PATH "$ENV{PNETCDF_PATH}")



