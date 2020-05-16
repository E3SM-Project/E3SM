# Allow a few variants for specifying the kokkos installation directory
if (Kokkos_DIR)
  SET (KOKKOS_INSTALL_DIR ${Kokkos_DIR})
elseif (KOKKOS_DIR)
  SET (KOKKOS_INSTALL_DIR ${KOKKOS_DIR})
else ()
  # Build kokkos submodule if user did not specify Kokkos_DIR.
  set(KOKKOS_SRC    ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/kokkos)
  set(KOKKOS_BINARY ${CMAKE_CURRENT_BINARY_DIR}/externals/kokkos)

  # We want Kokkos to be in debug mode if scream is in debug mode. This is a bit hacky
  # since we can't use KOKKOS_GMAKE_DEVICES yet.
  set(NO_DEBUG_ARGS Volta70 Pascal60)
  if (CMAKE_BUILD_TYPE_ci STREQUAL "debug" AND NOT KOKKOS_ARCH IN_LIST NO_DEBUG_ARGS)
    set(KOKKOS_ENABLE_DEBUG TRUE)
  endif()

  add_subdirectory(${KOKKOS_SRC} ${KOKKOS_BINARY})

  set (Kokkos_INCLUDE_DIR "${KOKKOS_SRC}/core/src;${KOKKOS_SRC}/algorithms/src;${KOKKOS_BINARY}")
  list(APPEND KOKKOS_LDFLAGS "-L${KOKKOS_BINARY}")
  list(REMOVE_ITEM KOKKOS_LDFLAGS "-lkokkos")
  list(APPEND LDFLAGS ${KOKKOS_LDFLAGS})
  set(KOKKOS_LIBS kokkos)
endif ()

if (KOKKOS_INSTALL_DIR)
  # Set some kokkos variables
  include (${KOKKOS_INSTALL_DIR}/kokkos_generated_settings.cmake)
  set (Kokkos_INCLUDE_DIR ${KOKKOS_INSTALL_DIR}/include)
  set (Kokkos_LIBRARY_DIR ${KOKKOS_INSTALL_DIR}/lib)
  set (Kokkos_LINK_FLAGS ${KOKKOS_LINK_FLAGS})

  # LB: is this really needed?
  string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
  set(KOKKOS_LDFLAGS_STR "")
  foreach(LDFLAG ${KOKKOS_LDFLAGS})
    if (${LDFLAG} IN_LIST KOKKOS_CXXFLAGS)
      message("Skipping repeat flag ${LDFLAG}")
    else()
      set(KOKKOS_LDFLAGS_STR "${KOKKOS_LDFLAGS_STR} ${LDFLAG}")
    endif()
  endforeach()

  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${KOKKOS_CXXFLAGS_STR}")
endif()
