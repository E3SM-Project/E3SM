include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

message(STATUS "tuo PROJECT_NAME=${PROJECT_NAME}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
	option(Kokkos_ARCH_ZEN4 "" ON)
	include(${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
endif()

# set(CMAKE_EXE_LINKER_FLAGS " -lxpmem  -L/opt/cray/pe/mpich/9.0.1/gtl/lib -lmpi_gtl_hsa -Wl,-rpath,/opt/cray/pe/mpich/9.0.1/gtl/lib") 
# set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

option(SCREAM_SMALL_KERNELS "Use small, non-monolothic kokkos kernels for ALL components that support them" ON)

message(STATUS "tuo CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
			set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
		endif()
	endif()
else()
	set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()
