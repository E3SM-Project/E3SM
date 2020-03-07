set (KOKKOS_INSTALLATION_NEEDED FALSE CACHE LOGICAL "")
if (NOT DEFINED E3SM_KOKKOS_PATH)
  # Build kokkos submodule if user did not specify KOKKOS_PATH.
  set (KOKKOS_SRC ${CMAKE_SOURCE_DIR}/../../externals/kokkos)
  set (KOKKOS_BUILD_DIR ${CMAKE_BINARY_DIR}/kokkos/build)

  # kokkos-containers is likely not needed, but it's built anyways,
  # so we may as well add its include path. kokkos-algorithms may
  # be used for the Kokkos_Random.hpp header.
  SET (KOKKOS_INCLUDE_DIR
       ${KOKKOS_SRC}/core/src
       ${KOKKOS_SRC}/containers/src
       ${KOKKOS_SRC}/algorithms/src
       ${KOKKOS_BUILD_DIR})
  SET (KOKKOS_LIBRARY_DIR ${KOKKOS_BUILD_DIR})
  message ("Using Kokkos built internally to E3SM, in ${KOKKOS_BUILD_DIR}")
  if (NOT E3SM_INTERNAL_KOKKOS_ALREADY_BUILT)
    # Nobody else already added kokkos subdirectory, so we will do it
    set (KOKKOS_INSTALLATION_NEEDED TRUE)
    set (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)
  endif()

else ()
  message ("The installation of kokkos in ${E3SM_KOKKOS_PATH} will be used.")
  SET (KOKKOS_INCLUDE_DIR ${E3SM_KOKKOS_PATH}/include)
  SET (KOKKOS_LIBRARY_DIR ${E3SM_KOKKOS_PATH}/lib64)
endif ()

SET (KOKKOS_LIBRARIES "kokkoscore;kokkoscontainers")
SET (KOKKOS_TPL_LIBRARIES "dl")

MESSAGE (STATUS "Kokkos installation found! Here are the details:")
MESSAGE ("    KOKKOS_INCLUDE_DIR: ${KOKKOS_INCLUDE_DIR}")
MESSAGE ("    KOKKOS_LIBRARY_DIR: ${KOKKOS_LIBRARY_DIR}")
MESSAGE ("    KOKKOS_LIBRARIES:   ${KOKKOS_LIBRARIES}")

macro(install_kokkos_if_needed)
  if (KOKKOS_INSTALLATION_NEEDED)
    if (ENABLE_OPENMP)
      set (Kokkos_ENABLE_OPENMP TRUE)
    endif ()
    add_subdirectory (${KOKKOS_SRC} ${KOKKOS_LIBRARY_DIR})
  endif ()
endmacro()

macro(link_to_kokkos targetName)
  # INCLUDE_DIRECTORIES("${KOKKOS_INCLUDE_DIR}")
  TARGET_INCLUDE_DIRECTORIES(${targetName} SYSTEM PUBLIC ${KOKKOS_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(${targetName} ${KOKKOS_TPL_LIBRARIES} ${KOKKOS_LIBRARIES} -L${KOKKOS_LIBRARY_DIR} ${LINKER_ADD_ON})
endmacro(link_to_kokkos)
