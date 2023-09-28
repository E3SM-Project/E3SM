# - Try to find Netcdf
#
# Placeholder code until spio gets a proper config.cmake
#
# Once done, this will define:
#
#   The "netcdf" target
#

# Pnetcdf is optional, and only if not running serial
if (NOT MPILIB STREQUAL mpi-serial)
  if (PNETCDF_PATH)
    find_library(pnetcdf_lib pnetcdf REQUIRED PATHS ${PNETCDF_PATH}/lib)
    find_path (pnetcdf_incdir pnetcdf.h REQUIRED PATHS ${PNETCDF_PATH}/include)
  endif()
endif()

if (NETCDF_C_PATH)
  # Sanity checks
  if (NOT NETCDF_FORTRAN_PATH)
    message(FATAL_ERROR "NETCDF_C_PATH specified without NETCDF_FORTRAN_PATH")
  endif()
  if (NOT EXISTS ${NETCDF_C_PATH}/lib AND NOT EXISTS ${NETCDF_C_PATH}/lib64)
    message(FATAL_ERROR "NETCDF_C_PATH does not contain a lib or lib64 directory")
  endif ()
  if (NOT EXISTS ${NETCDF_FORTRAN_PATH}/lib AND NOT EXISTS ${NETCDF_FORTRAN_PATH}/lib64)
    message(FATAL_ERROR "NETCDF_FORTRAN_PATH does not contain a lib or lib64 directory")
  endif ()

  # Find the libraries
  find_library(netcdf_c_lib netcdf  REQUIRED PATHS ${NETCDF_C_PATH}/lib ${NETCDF_C_PATH}/lib64)
  find_library(netcdf_f_lib netcdff REQUIRED PATHS ${NETCDF_FORTRAN_PATH}/lib ${NETCDF_FORTRAN_PATH}/lib64)
  find_path (netcdf_c_incdir netcdf.h REQUIRED PATHS ${NETCDF_C_PATH}/include)
  find_path (netcdf_f_incdir netcdf.inc REQUIRED PATHS ${NETCDF_FORTRAN_PATH}/include)

elseif (NETCDF_FORTRAN_PATH)
  message(FATAL_ERROR "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")
elseif (NETCDF_PATH)

  # Sanity checks
  if (NOT EXISTS ${NETCDF_PATH}/lib AND NOT EXISTS ${NETCDF_PATH}/lib64)
    message(FATAL_ERROR "NETCDF_PATH does not contain a lib or lib64 directory")
  endif ()

  find_library(netcdf_c_lib netcdf  REQUIRED PATHS ${NETCDF_PATH}/lib ${NETCDF_PATH}/lib64)
  find_library(netcdf_f_lib netcdff REQUIRED PATHS ${NETCDF_PATH}/lib ${NETCDF_PATH}/lib64)
  find_path (netcdf_c_incdir netcdf.h REQUIRED PATHS ${NETCDF_PATH}/include)
  find_path (netcdf_f_incdir netcdf.inc REQUIRED PATHS ${NETCDF_PATH}/include)
else()
  message(FATAL_ERROR "NETCDF not found: Define NETCDF_PATH or NETCDF_C_PATH and NETCDF_FORTRAN_PATH in config_machines.xml or config_compilers.xml")
endif()

set (pnetcdf_lib ${pnetcdf_lib})
set (netcdf_c_lib ${netcdf_c_lib})
set (netcdf_f_lib ${netcdf_f_lib})
set (pnetcdf_incdir ${pnetcdf_incdir})
set (netcdf_c_incdir ${netcdf_c_incdir})
set (netcdf_f_incdir ${netcdf_f_incdir})

# Create the interface library, and set target properties
add_library(netcdf INTERFACE)
target_link_libraries(netcdf INTERFACE ${pnetcdf_lib};${netcdf_c_lib};${netcdf_f_lib})
target_include_directories(netcdf INTERFACE ${pnetcdf_incdir};${netcdf_c_incdir};${netcdf_f_incdir})
if (pnetcdf_lib)
  target_compile_definitions(netcdf INTERFACE "-D_PNETCDF")
endif()
