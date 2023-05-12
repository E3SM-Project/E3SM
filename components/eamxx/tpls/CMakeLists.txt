# Here, we process a few shared libraries shipped with e3sm.
# For each of them, we either
#   - build the library (standalone build) or
#   - wrap pre-built library in a CMake target (CIME build)

# First pioc/piof, since we link against it in csm_share (at least in CIME build)
include (${SCREAM_BASE_DIR}/cmake/tpls/Scorpio.cmake)
CreateScorpioTargets()

# Then csm_share
include (${SCREAM_BASE_DIR}/cmake/tpls/CsmShare.cmake)
CreateCsmShareTarget()

if (SCREAM_CIME_BUILD)
  # For CIME runs, wrap mct in a target too
  include (${SCREAM_BASE_DIR}/cmake/tpls/Mct.cmake)
  CreateMctTarget()
endif()

# MAM aerosol support
if (SCREAM_ENABLE_MAM)
  # We use CMake's ExternalProject capability to build and install Haero.
  include(ExternalProject)

  # Build and install the Haero aerosol package interface.
  if (SCREAM_DOUBLE_PRECISION)
    set(HAERO_PRECISION "double")
  else()
    set(HAERO_PRECISION "single")
  endif()
  set(HAERO_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/externals/haero")
  set(HAERO_CMAKE_OPTS
    -DCMAKE_INSTALL_PREFIX=${HAERO_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DHAERO_ENABLE_GPU=${Kokkos_ENABLE_CUDA}
    -DHAERO_ENABLE_MPI=ON
    -DHAERO_PRECISION=${HAERO_PRECISION}
    -DEKAT_SOURCE_DIR=${PROJECT_SOURCE_DIR}/../../externals/ekat
    -DEKAT_BINARY_DIR=${PROJECT_BINARY_DIR}/externals/ekat)
  ExternalProject_Add(haero_proj
                      PREFIX ${PROJECT_BINARY_DIR}/externals/haero
                      SOURCE_DIR ${PROJECT_SOURCE_DIR}/../../externals/haero
                      BINARY_DIR ${PROJECT_BINARY_DIR}/externals/haero
                      CMAKE_ARGS ${HAERO_CMAKE_OPTS}
                      DEPENDS ekat ekat_test_session ekat_test_main
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)

  # Build and install MAM4xx.
  set(MAM4XX_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/externals/mam4xx")
  set(MAM4XX_CMAKE_OPTS
    -DCMAKE_INSTALL_PREFIX=${MAM4XX_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DMAM4XX_HAERO_DIR=${HAERO_INSTALL_PREFIX})
  ExternalProject_Add(mam4xx_proj
                      PREFIX ${PROJECT_BINARY_DIR}/externals/mam4xx
                      SOURCE_DIR ../../../../externals/mam4xx
                      BINARY_DIR ${PROJECT_BINARY_DIR}/externals/mam4xx
                      CMAKE_ARGS ${MAM4XX_CMAKE_OPTS}
                      DEPENDS haero_proj
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)

  # Bring in MAM4xx-related targets by including the generated mam4xx.cmake
  #  list(APPEND CMAKE_MODULE_PATH ${MAM4XX_INSTALL_PREFIX}/share)
  #  include(mam4xx)

endif()
