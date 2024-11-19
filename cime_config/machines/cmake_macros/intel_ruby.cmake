string(APPEND CPPDEFS " -DNO_SHR_VMATH -DCNL")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check all -ftrapuv")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/")

list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-cascadelake/intel-2021.6.0/hdf5-1.10.7-ewjpbjdhjgjzrzjcvwyjyuulaesbsjhg/lib")
list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-cascadelake/intel-2021.6.0/netcdf-c-4.4.1.1-vaxofekwvnvngh7wptmzkwdb7tkzvesn/lib")
list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-cascadelake/intel-2021.6.0/netcdf-fortran-4.4.4-3pzbx2unddhladhubaahhhysjmprzqi2/lib")
list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-cascadelake/intel-2021.6.0/parallel-netcdf-1.11.0-tzgdalakmem7tod6cruhqyeackeix5q5/lib")

set(KOKKOS_OPTIONS "--with-serial --ldflags='-L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/'")
