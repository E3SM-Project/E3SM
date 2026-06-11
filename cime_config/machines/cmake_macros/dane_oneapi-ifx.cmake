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

set(NETCDF_C_PATH "$ENV{NETCDF_C_PATH}")
set(NETCDF_FORTRAN_PATH "$ENV{NETCDF_FORTRAN_PATH}")
set(MKL_PATH "/usr/tce/packages/mkl/mkl-2022.1.0/lib/intel64")

list(APPEND CMAKE_BUILD_RPATH
	"${NETCDF_C_PATH}/lib"
	"${NETCDF_FORTRAN_PATH}/lib"
	"${MKL_PATH}"
)

list(APPEND CMAKE_INSTALL_RPATH
	"${NETCDF_C_PATH}/lib"
	"${NETCDF_FORTRAN_PATH}/lib"
	"${MKL_PATH}"
)

string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,--enable-new-dtags")
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--enable-new-dtags")

# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#     message("${_variableName}=${${_variableName}}")
# endforeach()
