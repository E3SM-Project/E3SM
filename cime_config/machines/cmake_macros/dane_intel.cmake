set(MPICC "/usr/tce/packages/mvapich2-tce/mvapich2-2.3.7-intel-2025.2.0/bin/mpicc")
set(MPICXX "/usr/tce/packages/mvapich2-tce/mvapich2-2.3.7-intel-2025.2.0/bin/mpicxx")
set(MPIFC "/usr/tce/packages/mvapich2-tce/mvapich2-2.3.7-intel-2025.2.0/bin/mpifort")
set(SCC "/usr/tce/packages/intel/intel-2025.2.0-magic/bin/icx")
set(SCXX "/usr/tce/packages/intel/intel-2025.2.0-magic/bin/icpx")
set(SFC "/usr/tce/packages/intel/intel-2025.2.0-magic/bin/ifx")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
	if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "2025.3.0")
		string(APPEND CMAKE_Fortran_FLAGS " -mllvm -disable-hir-temp-cleanup")
	endif()
endif()

set(MKL_PATH "/usr/tce/packages/mkl/mkl-2022.1.0/lib/intel64")

# Keep the runtime search path attached to the binaries instead of relying on
# LD_LIBRARY_PATH in the user's shell or batch environment.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

foreach (_prefix
	"$ENV{NETCDF_C_PATH}"
	"$ENV{NETCDF_FORTRAN_PATH}"
	"$ENV{PNETCDF_PATH}"
	"$ENV{HDF5_ROOT}"
	"$ENV{TEMPESTREMAP_ROOT}"
	"$ENV{MOAB_ROOT}")
	list(APPEND CMAKE_BUILD_RPATH "${_prefix}/lib" "${_prefix}/lib64")
	list(APPEND CMAKE_INSTALL_RPATH "${_prefix}/lib" "${_prefix}/lib64")
endforeach()

list(APPEND CMAKE_BUILD_RPATH "${MKL_PATH}")
list(APPEND CMAKE_INSTALL_RPATH "${MKL_PATH}")

list(REMOVE_DUPLICATES CMAKE_BUILD_RPATH)
list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)

if (NOT "${MPILIB}" STREQUAL "mpi-serial")
	find_library(E3SM_PNETCDF_LINK_LIB pnetcdf REQUIRED HINTS "$ENV{PNETCDF_PATH}/lib" NO_DEFAULT_PATH)
	# CMake already discovers PnetCDF through the netcdf interface target, but on this
	# platform the static Scorpio archive still needs libpnetcdf repeated late on the
	# final Fortran link line so the linker can resolve ncmpi_* symbols from libpioc.a.
	string(APPEND CMAKE_Fortran_STANDARD_LIBRARIES " ${E3SM_PNETCDF_LINK_LIB}")
endif()

string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,--enable-new-dtags")
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--enable-new-dtags")

# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#     message("${_variableName}=${${_variableName}}")
# endforeach()
