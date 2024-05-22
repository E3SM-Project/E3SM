include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
# kokkos sycl is on in the above file
#include (${EKAT_MACH_FILES_PATH}/kokkos/sycl.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

#AB flags from ekat
# -fsycl -fsycl-unnamed-lambda -sycl-std=2020 -qopenmp-simd -Wsycl-strict -fsycl-device-code-split=per_kernel

#SET(MPICH_DIR "/soft/restricted/CNDA/updates/mpich/52.2/mpich-ofi-all-icc-default-pmix-gpu-drop52/" CACHE STRING "")

set(CMAKE_CXX_FLAGS " -\-intel -Xclang -fsycl-allow-virtual-functions -mlong-double-64 -O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "-fc=ifx -O3 -DNDEBUG -DCPRINTEL -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS " -lifcore -\-intel -Xclang -mlong-double-64 -DNDEBUG  -fortlib" CACHE STRING "" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS " -Wl,-\-defsym,main=MAIN_\_ -lifcore -\-intel -Xclang -fsycl-allow-virtual-functions -lsycl -mlong-double-64 -DNDEBUG ${SYCL_LINK_FLAGS} -fortlib  -L${MPICH_DIR}/lib" CACHE STRING "" FORCE)



#
#        <env name="NETCDF_C_PATH">/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2023.05.15.007</env>
#        <env name="NETCDF_FORTRAN_PATH">/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2023.05.15.007</env>
#        <env name="PNETCDF_PATH">/lus/gecko/projects/CSC249ADSE15_CNDA/software/pnetcdf/1.12.3/oneapi.eng.2023.05.15.007</env>


#this is needed for cime builds!
set(NETCDF_PATH "/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2023.05.15.007")
set(NETCDF_DIR  "/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2023.05.15.007")
set(NETCDF_C_PATH "/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2023.05.15.007")
set(NETCDF_C      "/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2023.05.15.007")
#this one is for rrtmgp
set(NetCDF_C_PATH "/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-c/4.9.2/oneapi.eng.2023.05.15.007" CACHE STRING "")
set(NETCDF_FORTRAN_PATH "/lus/gecko/projects/CSC249ADSE15_CNDA/software/netcdf-fortran/4.6.1/oneapi.eng.2023.05.15.007")
set(PNETCDF_PATH "/lus/gecko/projects/CSC249ADSE15_CNDA/software/pnetcdf/1.12.3/oneapi.eng.2023.05.15.007")



