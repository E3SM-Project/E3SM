string(APPEND CPPDEFS " -DNO_SHR_VMATH -DCNL")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check all -ftrapuv")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/")

list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-sapphirerapids/intel-2021.6.0/hdf5-1.10.7-766kapalbrdntu2pcgdgbhg2ch26gsuv/lib")
list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-sapphirerapids/intel-2021.6.0/netcdf-c-4.4.1.1-2uznnlwgiezxute6iyqzqjrpolokeaib/lib")
list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-sapphirerapids/intel-2021.6.0/netcdf-fortran-4.4.4-itpstyordbern7vlulmlnt47eeeokzfp/lib")
list(APPEND CMAKE_BUILD_RPATH "/usr/workspace/e3sm/spack/libs/linux-rhel8-sapphirerapids/intel-2021.6.0/parallel-netcdf-1.11.0-26sxm4mormsglmhi24poix7sugbigkck/lib")

set(KOKKOS_OPTIONS "--with-serial --ldflags='-L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/'")
