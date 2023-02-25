set (KOKKOS_LIBRARIES "kokkoscore;kokkoscontainers")
set (KOKKOS_TPL_LIBRARIES "dl")

set (KOKKOS_INSTALLATION_NEEDED FALSE CACHE BOOL "")
if (DEFINED E3SM_KOKKOS_PATH)
  message (STATUS "The E3SM installation of kokkos in ${E3SM_KOKKOS_PATH} will be used. Here are the details:")
  set (KOKKOS_INCLUDE_DIR ${E3SM_KOKKOS_PATH}/include)
  set (KOKKOS_LIBRARY_DIR ${E3SM_KOKKOS_PATH}/lib64)

  message ("    KOKKOS_INCLUDE_DIR: ${KOKKOS_INCLUDE_DIR}")
  message ("    KOKKOS_LIBRARY_DIR: ${KOKKOS_LIBRARY_DIR}")
  message ("    KOKKOS_LIBRARIES:   ${KOKKOS_LIBRARIES}")
  set (KOKKOS_LIBS_ARE_TARGETS FALSE)
elseif (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT)
  message (STATUS "All kokkos libraries needed by Homme are already built by this cmake project.\n")
  set (KOKKOS_LIBS_ARE_TARGETS TRUE)
else()
  # As last resort, build kokkos submodule
  set (KOKKOS_SRC ${CMAKE_SOURCE_DIR}/../../externals/kokkos)
  set (KOKKOS_BUILD_DIR ${CMAKE_BINARY_DIR}/kokkos/build)

  # kokkos-containers is likely not needed, but it's built anyways,
  # so we may as well add its include path. kokkos-algorithms may
  # be used for the Kokkos_Random.hpp header.
  set (KOKKOS_INCLUDE_DIR
       ${KOKKOS_SRC}/core/src
       ${KOKKOS_SRC}/containers/src
       ${KOKKOS_SRC}/algorithms/src
       ${KOKKOS_BUILD_DIR})
  set (KOKKOS_LIBRARY_DIR ${KOKKOS_BUILD_DIR})
  message ("Building Kokkos (from E3SM internal submodule) in ${KOKKOS_BUILD_DIR}")

  # Make a note that we need to add the kokkos subdir (in the cmake sense)
  set (KOKKOS_INSTALLATION_NEEDED TRUE)

  # Notify that E3SM's internal kokkos is built
  set (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)

  # Make a note that, since kokkos is built internally, it's a valid cmake target
  # (in particular, we don't need to set compile/include/link flags)
  set (KOKKOS_LIBS_ARE_TARGETS TRUE)
endif ()

macro(install_kokkos_if_needed)
  if (KOKKOS_INSTALLATION_NEEDED)
    if (ENABLE_OPENMP)
      set (Kokkos_ENABLE_OPENMP TRUE)
    endif ()
    add_subdirectory (${KOKKOS_SRC} ${KOKKOS_LIBRARY_DIR})
  endif ()
endmacro()

macro(link_to_kokkos targetName)
  target_include_directories(${targetName} SYSTEM PUBLIC ${KOKKOS_INCLUDE_DIR})
  if (KOKKOS_LIBS_ARE_TARGETS)
    target_link_libraries (${targetName} ${KOKKOS_LIBRARIES})
  else()
    target_link_libraries(${targetName} ${KOKKOS_TPL_LIBRARIES} ${KOKKOS_LIBRARIES} -L${KOKKOS_LIBRARY_DIR})
  endif()
endmacro(link_to_kokkos)
