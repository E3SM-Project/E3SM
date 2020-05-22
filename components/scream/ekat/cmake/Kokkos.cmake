# Define a global property to check if Kokkos has already been built
define_property(GLOBAL
                PROPERTY EKAT_KOKKOS_BUILT
                BRIEF_DOCS "Wheter kokkos subdir has already been processed"
                FULL_DOCS "This property is used by cmake to ensure that Kokkos
                           submodule directory is only processed once (with add_subdirectory).")

get_property(IS_EKAT_KOKKOS_BUILT GLOBAL PROPERTY EKAT_KOKKOS_BUILT SET)

# Process the kokkos source directory
if (NOT IS_EKAT_KOKKOS_BUILT)

  if (NOT Kokkos_SOURCE_DIR)
    message (FATAL_ERROR "Error! Please, specify path to Kokkos in Kokkos_SOURCE_DIR.\n")
  elseif (NOT EXISTS ${Kokkos_SOURCE_DIR})
    message (FATAL_ERROR "Error! Please specify a valid source folder for kokkos.\n"
                         "       Provided path: ${Kokkos_SOURCE_DIR}")
  endif()
  set(Kokkos_BINARY_DIR ${CMAKE_BINARY_DIR}/externals/kokkos)

  # We want Kokkos to be in debug mode if the host project is in debug mode. This is a bit hacky
  # since we can't use Kokkos_GMAKE_DEVICES yet.
  set(NO_DEBUG_ARGS Volta70 Pascal60)
  if (CMAKE_BUILD_TYPE_ci STREQUAL "debug" AND NOT Kokkos_ARCH IN_LIST NO_DEBUG_ARGS)
    set(Kokkos_ENABLE_DEBUG TRUE)
  endif()

  add_subdirectory(${Kokkos_SOURCE_DIR} ${Kokkos_BINARY_DIR})

  set (Kokkos_INCLUDE_DIRS
    ${Kokkos_SOURCE_DIR}/core/src
    ${Kokkos_SOURCE_DIR}/algorithms/src
    ${Kokkos_BINARY_DIR}
  )
  set(Kokkos_LIBRARIES kokkos)

  # Make sure it is processed only once
  set_property(GLOBAL PROPERTY EKAT_KOKKOS_BUILT TRUE)
endif ()
