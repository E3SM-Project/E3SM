# This should be included by all machine files and common_setup
# should be called immediately after.
#
# The macro sets EKAT_MACH_FILES_PATH which can then be used by
# individual machine files to include EKAT machine/kokkos files.
# It also sets the critical SCREAM_MACHINE variable.
#
# The SCREAM_MACHINE does not have to be known to EKAT, but more
# work is needed in scream's $machine.cmake file if it is not
# (the kokkos and mpi settings will need to be manually set).
macro(common_setup)
  set(EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

  # Set SCREAM_MACHINE if it's not set.
  if (NOT SCREAM_MACHINE)
    if (DEFINED ENV{SCREAM_MACHINE})
      set(SCREAM_MACHINE $ENV{SCREAM_MACHINE} CACHE STRING "")
    elseif (MACH)
      set(SCREAM_MACHINE ${MACH} CACHE STRING "")
    else()
      get_filename_component(TEMP ${CMAKE_CURRENT_LIST_FILE} NAME_WLE)
      set(SCREAM_MACHINE ${TEMP} CACHE STRING "")
    endif()
  endif()

  # include EKAT machine file if it exists
  set(EKAT_MACH_FILE_PATH ${EKAT_MACH_FILES_PATH}/${SCREAM_MACHINE}.cmake)
  if (EXISTS ${EKAT_MACH_FILE_PATH})
    include(${EKAT_MACH_FILE_PATH})
  endif()

endmacro()
