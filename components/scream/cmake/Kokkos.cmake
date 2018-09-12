# Allow a few variants for specifying the kokkos installation directory
if (Kokkos_DIR)
  SET (KOKKOS_INSTALL_DIR ${Kokkos_DIR})
elseif (KOKKOS_DIR)
  SET (KOKKOS_INSTALL_DIR ${KOKKOS_DIR})
else ()
  message (FATAL_ERROR "SCREAM requires Kokkos. Set one of Kokkos_DIR or KOKKOS_DIR to point to the base directory of the Kokkos installation.")
endif ()

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
