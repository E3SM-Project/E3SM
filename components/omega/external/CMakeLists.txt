# Add external packages

# Add the spdlog library
add_subdirectory(
  ${E3SM_EXTERNALS_ROOT}/ekat/extern/spdlog
  ${CMAKE_CURRENT_BINARY_DIR}/ekat/extern/spdlog
)

# Add the YAKL library
add_subdirectory(
  ${E3SM_EXTERNALS_ROOT}/YAKL
  ${CMAKE_CURRENT_BINARY_DIR}/YAKL
)

if (CMAKE_VERSION VERSION_GREATER "3.18.0")
  set_property(TARGET yakl PROPERTY CUDA_ARCHITECTURES OFF)
endif()

if(OMEGA_DEBUG)
  target_compile_definitions(yakl INTERFACE YAKL_DEBUG=1)
endif()

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

# Add the Scorpio library
add_subdirectory(
  ${E3SM_EXTERNALS_ROOT}/scorpio
  ${CMAKE_CURRENT_BINARY_DIR}/scorpio
)

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
