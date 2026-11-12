# Add external packages

include(ExternalProject)

# Add GSW-C project
ExternalProject_Add(
    gswteos
    BINARY_DIR ${OMEGA_BUILD_DIR}/external/GSW-C
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory
		${OMEGA_SOURCE_DIR}/external/GSW-C
		${OMEGA_BUILD_DIR}/external/GSW-C
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} all
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${OMEGA_BUILD_DIR}/external/GSW-C/libgswteos-10.so
		${OMEGA_BUILD_DIR}/external/GSW-C/gsw_check
)

add_library(gswteos-10 SHARED IMPORTED)
set_target_properties(gswteos-10 PROPERTIES
    IMPORTED_LOCATION ${OMEGA_BUILD_DIR}/external/GSW-C/libgswteos-10.so
)

# Add the spdlog library
# Guard on both of bare and namespaced names.
if (NOT TARGET spdlog AND NOT TARGET spdlog::spdlog)
	add_subdirectory(
	  ${OMEGA_SOURCE_DIR}/external/spdlog
	  ${CMAKE_CURRENT_BINARY_DIR}/external/spdlog
	)
        # workaround for craygnu-mphipcc error on frontier when compiling spdlog
        if("${OMEGA_ARCH}" STREQUAL "HIP")
          target_compile_definitions(spdlog PRIVATE NDEBUG=1)
        endif()
endif()

# Omega links the bare name (src/CMakeLists.txt), so guarantee it exists.
if (NOT TARGET spdlog)
	add_library(spdlog INTERFACE)
	target_link_libraries(spdlog INTERFACE spdlog::spdlog)
endif()

# Add the yaml-cpp library
if (NOT TARGET yaml-cpp AND NOT TARGET yaml-cpp::yaml-cpp)
	add_subdirectory(
	  ${OMEGA_SOURCE_DIR}/external/yaml-cpp
	  ${CMAKE_CURRENT_BINARY_DIR}/external/yaml-cpp
	)
endif()

# Omega links the bare name (src/CMakeLists.txt)
if (NOT TARGET yaml-cpp)
	add_library(yaml-cpp INTERFACE)
	target_link_libraries(yaml-cpp INTERFACE yaml-cpp::yaml-cpp)
endif()

# Add the Kokkos library.
# Reuse-if-present: when Omega is built alongside EAMxx (or a shared EKAT),
# the Kokkos::kokkos target already exists with the correct per-machine
# architecture and compiler configuration, so Omega simply reuses it and this
# block is skipped entirely. Otherwise Omega must build Kokkos itself.
if (NOT TARGET Kokkos::kokkos AND NOT TARGET kokkos)
	if ("${OMEGA_BUILD_MODE}" STREQUAL "E3SM")
		list(APPEND CMAKE_MODULE_PATH
			${E3SM_EXTERNALS_ROOT}/ekat/cmake
			${E3SM_EXTERNALS_ROOT}/ekat/cmake/tpls
		)
		include(EkatBuildKokkos)
	else()
		add_subdirectory(
			${E3SM_EXTERNALS_ROOT}/ekat/extern/kokkos
			${CMAKE_CURRENT_BINARY_DIR}/ekat/extern/kokkos
		)
	endif()
endif()

# Omega links the namespaced name (src/CMakeLists.txt).
if (NOT TARGET Kokkos::kokkos)
	add_library(Kokkos::kokkos ALIAS kokkos)
endif()

# Link the case's own Scorpio and GPTL in a coupled build.
# E3SM has already built both into the case's shared libraries, and
# components/cmake/modules/FindPIO.cmake wraps them in one interface target,
# "spio", with no target named "pioc" or "gptl". Without the two wrappers
# below, the guards that follow find nothing to reuse: Omega builds a second
# Scorpio, and with it Scorpio's bundled GPTL, so e3sm.exe links two copies of
# each and takes GPTL from the private one, which lacks -DHAVE_SLASHPROC.
# EAMxx wraps "spio" the same way (components/eamxx/cmake/tpls/Scorpio.cmake);
# when it has, these skip.
if (TARGET spio AND NOT TARGET pioc)
	add_library(pioc INTERFACE)
	target_link_libraries(pioc INTERFACE spio)
	# pio.h is in the source tree; the generated pio_config.h it includes is
	# only in the shared include directory.
	target_include_directories(pioc INTERFACE
	  ${E3SM_EXTERNALS_ROOT}/scorpio/src/clib
	  ${INSTALL_SHAREDPATH}/include
	)
endif()

if (TARGET spio AND NOT TARGET gptl)
	add_library(gptl INTERFACE)
	target_link_libraries(gptl INTERFACE spio)
	target_include_directories(gptl INTERFACE ${INSTALL_SHAREDPATH}/include)
	# gptl.h declares its MPI_Comm signatures only under HAVE_MPI, which
	# share/timing builds the case's GPTL with unless MPILIB is mpi-serial.
	# Scorpio's bundled GPTL exported it PUBLIC, and pacer relies on it:
	# without it, passing an MPI_Comm fails to compile wherever MPI_Comm is
	# not an int, as with OpenMPI.
	if (NOT MPILIB STREQUAL "mpi-serial")
		target_compile_definitions(gptl INTERFACE HAVE_MPI)
	endif()
endif()

# Add the Scorpio library
# Deliberately left as a single-name guard: unlike spdlog and yaml-cpp, "pioc"
# has no namespaced alias to also test, so one name already covers every path
# that supplies it.
if (NOT TARGET pioc)

	# forward env. variables to Scorpio build
	if(DEFINED ENV{NETCDF_PATH})
	  set(NetCDF_PATH $ENV{NETCDF_PATH}
			CACHE STRING "Path to NETCDF library" FORCE)
	endif()
	if(DEFINED ENV{NETCDF_C_PATH})
	  set(NetCDF_C_PATH $ENV{NETCDF_C_PATH}
			CACHE STRING "Path to NETCDF-C library" FORCE)
	endif()
	if(DEFINED ENV{NETCDF_FORTRAN_PATH})
	  set(NetCDF_Fortran_PATH $ENV{NETCDF_FORTRAN_PATH}
			CACHE STRING "Path to NETCDF-Fortran library" FORCE)
	endif()
	if(DEFINED ENV{PNETCDF_PATH})
	  set(PnetCDF_PATH $ENV{PNETCDF_PATH}
			CACHE STRING "Path to PNETCDF library" FORCE)
	endif()

	option(PIO_ENABLE_TOOLS "" OFF)

	# Locally shadow ADIOS2_FOUND, which leaks in from E3SM's top-level
	# find_package. Workaround for scorpio/src/clib/CMakeLists.txt:498, which
	# should test WITH_ADIOS2 (OFF here) instead of ADIOS2_FOUND.
	set(ADIOS2_FOUND FALSE)

	add_subdirectory(
	  ${E3SM_EXTERNALS_ROOT}/scorpio
	  ${CMAKE_CURRENT_BINARY_DIR}/scorpio
	)
endif()

# Add E3SM GPTL library
# Also a single-name guard: "gptl" has no namespaced alias anywhere.
if (NOT TARGET gptl)
	add_subdirectory(
	  ${E3SM_ROOT}/share/timing
	  ${CMAKE_CURRENT_BINARY_DIR}/gptl
	)
elseif (compile_threaded AND NOT TARGET spio)
	# The Scorpio add_subdirectory above already created a target named "gptl"
	# (externals/scorpio/src/gptl, pulled in whenever Scorpio's find_package(GPTL)
	# comes up empty), so the guard skipped E3SM's share/timing GPTL -- the only
	# one that sets -DTHREADED_OMP.
	#
	# Applies only when Omega built Scorpio itself, which in practice means a
	# standalone build. In a coupled build "gptl" is an interface wrapper
	# around "spio" -- the one above, or EAMxx's -- with nothing of its own
	# to patch, and the case's GPTL behind it is already threaded:
	# share/timing adds -DTHREADED_OMP whenever compile_threaded is set.
	#
	# Build Scorpio's copy threaded so it matches how E3SM uses it. Note the
	# Omega subtree does not inherit Macros.cmake's flags (CMAKE_C_FLAGS and
	# CMAKE_Fortran_FLAGS are empty here), so -fopenmp has to come along too.
	#
	# Ask for OpenMP's Fortran component only where Fortran is enabled in this
	# directory. A standalone build enables C and CXX only; Scorpio enables
	# Fortran in its own subdirectory, which this one does not see, so asking
	# for it here fails configure. Standalone Omega calls GPTL from C++ only.
	set(_gptl_omp_langs C)
	if (CMAKE_Fortran_COMPILER_LOADED)
		list(APPEND _gptl_omp_langs Fortran)
	endif()
	find_package(OpenMP REQUIRED COMPONENTS ${_gptl_omp_langs})
	target_compile_definitions(gptl PRIVATE THREADED_OMP)
	foreach(_lang IN LISTS _gptl_omp_langs)
		target_link_libraries(gptl PRIVATE OpenMP::OpenMP_${_lang})
	endforeach()
endif()

# Add the Pacer library
# Single-name guard: no E3SM-side "pacer" target exists anywhere, so Omega
# always builds it.
if (NOT TARGET pacer)
    add_library(pacer ${E3SM_ROOT}/share/pacer/Pacer.cpp)

    target_include_directories(
        pacer
        PRIVATE
        ${E3SM_ROOT}/share/timing
    )

    target_link_libraries(
        pacer
        PRIVATE
        gptl
        Kokkos::kokkos
    )

    target_compile_definitions(
        pacer
        PRIVATE
        PACER_HAVE_KOKKOS=1
     )
    if (OMEGA_EXTERNAL_PROF)
      target_compile_definitions(
          pacer
          PRIVATE
          PACER_ADD_RANGES=1
      )
    endif()
endif()

# Add the parmetis and related libraries

if(Parmetis_FOUND)
  add_library(parmetis STATIC IMPORTED GLOBAL)
  set_target_properties(parmetis PROPERTIES
     IMPORTED_LOCATION ${Parmetis_LIBRARY})
endif()

if(Metis_FOUND)
  add_library(metis STATIC IMPORTED GLOBAL)
  set_target_properties(metis PROPERTIES
     IMPORTED_LOCATION ${Metis_LIBRARY})
endif()

if(GKlib_FOUND)
  add_library(gklib STATIC IMPORTED GLOBAL)
  set_target_properties(gklib PROPERTIES
     IMPORTED_LOCATION ${GKlib_LIBRARY})
endif()

# Add the cpptrace library
# Note the real target is named "cpptrace-lib", not "cpptrace"
# (external/cpptrace/CMakeLists.txt sets target_name to avoid colliding with
# its own PROJECT_NAME); cpptrace::cpptrace is an ALIAS of it.
if (NOT TARGET cpptrace::cpptrace AND NOT TARGET cpptrace-lib)
	# Unset PACKAGE_VERSION before adding cpptrace to prevent ADIOS2's
	# PACKAGE_VERSION (e.g. 2.10.2) from leaking into libdwarf's configure_file
	unset(PACKAGE_VERSION)
	option(CPPTRACE_INHERIT_HOST_STANDARD "" ON)
        option(CPPTRACE_DISABLE_CXX_20_MODULES "" ON)
	add_subdirectory(
	  ${OMEGA_SOURCE_DIR}/external/cpptrace
	  ${CMAKE_CURRENT_BINARY_DIR}/external/cpptrace
	)
endif()

# Omega links the namespaced name (src/CMakeLists.txt).
if (NOT TARGET cpptrace::cpptrace)
	add_library(cpptrace::cpptrace ALIAS cpptrace-lib)
endif()
