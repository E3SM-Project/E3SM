# - Try to find Netcdf
#
# Placeholder code until spio gets a proper config.cmake
#
# Once done, this will define:
#
#   The "netcdf" target
#

if (TARGET netcdf)
  return()
endif()

# Sets GET_NETCDF_LIBS_RESULT
function(get_netcdf_libs ncpath nfpath)

  set(ncconfig ${ncpath}/bin/nc-config)
  set(nfconfig ${nfpath}/bin/nf-config)

  # Get C libs
  if (EXISTS ${ncconfig})
    execute_process(COMMAND ${ncconfig} --libs OUTPUT_VARIABLE nclibs OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()

  # Fall back to find_library
  if (NOT nclibs)
    find_library(nclibs_temp netcdf REQUIRED HINTS ${ncpath}/lib ${ncpath}/lib64)
    set(nclibs ${nclibs_temp})
  endif()

  # Get fortran libs
  if (EXISTS ${nfconfig})
    execute_process(COMMAND ${nfconfig} --flibs OUTPUT_VARIABLE nflibs OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()

  if (NOT nflibs)
    find_library(nflibs_temp netcdff REQUIRED HINTS ${nfpath}/lib ${nfpath}/lib64)
    set(nflibs ${nflibs_temp})
  endif()

  # C libs need to come last
  set(GET_NETCDF_LIBS_RESULT ${nflibs} ${nclibs} PARENT_SCOPE)
endfunction()

function(create_netcdf_target)

  # Grab things from env
  set(PNETCDF_PATH        $ENV{PNETCDF_PATH})
  set(NETCDF_PATH         $ENV{NETCDF_PATH})
  set(NETCDF_C_PATH       $ENV{NETCDF_C_PATH})
  set(NETCDF_FORTRAN_PATH $ENV{NETCDF_FORTRAN_PATH})

  # Pnetcdf is optional, and only if not running serial
  if (NOT MPILIB STREQUAL mpi-serial)
    if (PNETCDF_PATH)
      find_library(pnetcdf_lib pnetcdf REQUIRED HINTS ${PNETCDF_PATH}/lib)
      find_path (pnetcdf_incdir pnetcdf.h REQUIRED HINTS ${PNETCDF_PATH}/include)
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

    get_netcdf_libs(${NETCDF_C_PATH} ${NETCDF_FORTRAN_PATH})
    find_path (netcdf_c_incdir netcdf.h REQUIRED HINTS ${NETCDF_C_PATH}/include)
    find_path (netcdf_f_incdir netcdf.inc REQUIRED HINTS ${NETCDF_FORTRAN_PATH}/include)

  elseif (NETCDF_FORTRAN_PATH)
    message(FATAL_ERROR "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")

  elseif (NETCDF_PATH)
    # Sanity checks
    if (NOT EXISTS ${NETCDF_PATH}/lib AND NOT EXISTS ${NETCDF_PATH}/lib64)
      message(FATAL_ERROR "NETCDF_PATH does not contain a lib or lib64 directory")
    endif ()

    get_netcdf_libs(${NETCDF_PATH} ${NETCDF_PATH})
    find_path(netcdf_c_incdir netcdf.h REQUIRED HINTS ${NETCDF_PATH}/include)
    find_path(netcdf_f_incdir netcdf.inc REQUIRED HINTS ${NETCDF_PATH}/include)

  else()
    message(FATAL_ERROR "NETCDF not found: Define NETCDF_PATH or NETCDF_C_PATH and NETCDF_FORTRAN_PATH in config_machines.xml or config_compilers.xml")
  endif()

  set(pnetcdf_lib ${pnetcdf_lib})
  set(pnetcdf_incdir ${pnetcdf_incdir})
  set(netcdf_c_incdir ${netcdf_c_incdir})
  set(netcdf_f_incdir ${netcdf_f_incdir})

  # Create the interface library, and set target properties
  add_library(netcdf INTERFACE)
  target_link_libraries(netcdf INTERFACE ${pnetcdf_lib} ${GET_NETCDF_LIBS_RESULT})
  target_include_directories(netcdf INTERFACE ${pnetcdf_incdir};${netcdf_c_incdir};${netcdf_f_incdir})
endfunction()

create_netcdf_target()
