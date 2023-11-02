# If this is a CIME build, create IMPORTED target to wrap scorpio libs.
# Otherwise, simply add scorpio subdirectory.
set (E3SM_EXTERNALS_DIR ${CMAKE_CURRENT_LIST_DIR}/../../../../externals CACHE INTERNAL "")

set (SCREAM_TPLS_MODULE_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")

macro (CreateScorpioTargets)

  # Sanity check
  if (TARGET pioc OR TARGET piof)
    # We should not call this macro twice
    message (FATAL_ERROR "The Scorpio targets were already created!")
  endif()

  if (SCREAM_CIME_BUILD)
    find_package(PIO REQUIRED)

    add_library (pioc INTERFACE)
    target_link_libraries (pioc INTERFACE spio)

    add_library (piof INTERFACE)
    target_link_libraries (piof INTERFACE spio)

  else ()
    # Not a CIME build. We'll add scorpio as a subdir

    # We don't need (yet) SCORPIO tools
    option (PIO_ENABLE_TOOLS "Enable SCORPIO tools" OFF)

    # We want to use GPTL internally
    option (PIO_ENABLE_TIMING    "Enable the use of the GPTL timing library" ON)

    # This is the default, but just in case scorpio changes it
    option (PIO_ENABLE_FORTRAN "Enable the Fortran library builds" ON)

    add_subdirectory (${E3SM_EXTERNALS_DIR}/scorpio ${CMAKE_BINARY_DIR}/externals/scorpio)
    EkatDisableAllWarning(pioc)
    EkatDisableAllWarning(piof)
    EkatDisableAllWarning(gptl)
  endif ()
endmacro()
