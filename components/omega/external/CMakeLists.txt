# Add external packages

# Add the spdlog library
if (NOT TARGET spdlog::spdlog)
	add_subdirectory(
	  ${E3SM_EXTERNALS_ROOT}/ekat/extern/spdlog
	  ${CMAKE_CURRENT_BINARY_DIR}/ekat/extern/spdlog
	)
endif()

# Add the yaml-cpp library
if (NOT TARGET yaml-cpp::yaml-cpp)
	add_subdirectory(
	  ${E3SM_EXTERNALS_ROOT}/ekat/extern/yaml-cpp
	  ${CMAKE_CURRENT_BINARY_DIR}/ekat/extern/yaml-cpp
	)
endif()

# Add the Kokkos library
if (NOT TARGET Kokkos::kokkos)
	add_subdirectory(
      ${E3SM_EXTERNALS_ROOT}/ekat/extern/kokkos
      ${CMAKE_CURRENT_BINARY_DIR}/ekat/extern/kokkos
	)
endif()

# Add the Scorpio library
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

	add_subdirectory(
	  ${E3SM_EXTERNALS_ROOT}/scorpio
	  ${CMAKE_CURRENT_BINARY_DIR}/scorpio
	)
endif()

# Add E3SM GPTL library
if (NOT TARGET gptl)
	add_subdirectory(
	  ${E3SM_ROOT}/share/timing
	  ${CMAKE_CURRENT_BINARY_DIR}/gptl
	)
endif()

# Add the Pacer library
if (NOT TARGET pacer)
    add_library(pacer ${E3SM_ROOT}/share/pacer/Pacer.cpp)

    target_include_directories(
        pacer
        PRIVATE
        ${E3SM_ROOT}/share/timing
    )

    target_link_libraries(
        pacer
        PUBLIC
        gptl
    )
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
