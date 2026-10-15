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
    find_package(SCORPIO REQUIRED)
  else ()
    # Not a CIME build. We'll add scorpio as a subdir

    # We don't need (yet) SCORPIO tools
    option (PIO_ENABLE_TOOLS "Enable SCORPIO tools" OFF)

    # We want to use GPTL internally
    option (PIO_ENABLE_TIMING    "Enable the use of the GPTL timing library" ON)

    # This is the default, but just in case scorpio changes it
    option (PIO_ENABLE_FORTRAN "Enable the Fortran library builds" ON)

    # Use malloc (instead of the custom bget) for mallocs
    option (PIO_USE_MALLOC "Use native malloc (instead of bget package)" ON)

    add_subdirectory (${E3SM_EXTERNALS_DIR}/scorpio ${CMAKE_BINARY_DIR}/externals/scorpio)

    set (SCORPIO_Fortran_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/externals/scorpio/src/flib CACHE INTERNAL "SCORPIO Fortran include dirs")
    set (C_INCLUDE_DIRS "${E3SM_EXTERNALS_DIR}/scorpio/src/clib")
    list(APPEND C_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/externals/scorpio/src/clib")
    set (SCORPIO_C_INCLUDE_DIRS "${C_INCLUDE_DIRS}" CACHE INTERNAL "SCORPIO C include dirs")

    # Add include dirs to the GPTL target that SCORPIO builds (scorpio sets them as PRIVATE properties, unfortunately)
    if (NOT GPTL_PATH)
      set (GPTL_PATH ${E3SM_EXTERNALS_DIR}/scorpio/src/gptl CACHE INTERNAL "Path to GPTL library")
    endif ()
    target_include_directories(gptl PUBLIC $<BUILD_INTERFACE:${GPTL_PATH}>)
    target_include_directories(gptl PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/externals/scorpio/src/gptl>)

    EkatDisableAllWarning(pioc)
    EkatDisableAllWarning(piof)
    EkatDisableAllWarning(gptl)
  endif ()
endmacro()
