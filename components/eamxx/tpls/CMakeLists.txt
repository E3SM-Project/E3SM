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

if (NOT SCREAM_LIB_ONLY)
  include(BuildCprnc)
  BuildCprnc()
endif()

# MAM aerosol support
if (SCREAM_ENABLE_MAM)
  # We use CMake's ExternalProject capability to build and install Haero.
  include(ExternalProject)

  if (SCREAM_CIME_BUILD)
    # PROJECT_SOURCE_DIR is SCREAM_ROOT/components
    set(EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/../externals")
  else()
    # PROJECT_SOURCE_DIR is SCREAM_ROOT/components/eamxx
    set(EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/../../externals")
  endif()

  # Normalize CMAKE_BUILD_TYPE.
  string(TOLOWER "${CMAKE_BUILD_TYPE}" lc_build_type)
  if(${lc_build_type} STREQUAL "release")
    set(mam_build_type "Release")
  else()
    set(mam_build_type "Debug") # when in doubt...
  endif()

  # Build and install the Haero aerosol package interface.
  if (SCREAM_DOUBLE_PRECISION)
    set(HAERO_PRECISION "double")
  else()
    set(HAERO_PRECISION "single")
  endif()
  set(HAERO_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/externals/haero")
  set(HAERO_CMAKE_OPTS
    -DCMAKE_INSTALL_PREFIX=${HAERO_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${mam_build_type}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DHAERO_ENABLE_GPU=${Kokkos_ENABLE_CUDA}
    -DHAERO_ENABLE_MPI=ON
    -DHAERO_ENABLE_TESTS=OFF
    -DHAERO_PRECISION=${HAERO_PRECISION}
    -DEKAT_SOURCE_DIR=${EXTERNALS_DIR}/ekat
    -DEKAT_BINARY_DIR=${PROJECT_BINARY_DIR}/externals/ekat
    -DBUILD_SHARED_LIBS=OFF)
  ExternalProject_Add(haero_proj
                      PREFIX ${PROJECT_BINARY_DIR}/externals/haero
                      SOURCE_DIR ${EXTERNALS_DIR}/haero
                      BINARY_DIR ${PROJECT_BINARY_DIR}/externals/haero
                      CMAKE_ARGS ${HAERO_CMAKE_OPTS}
                      DEPENDS ekat ekat_test_session ekat_test_main
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE
                      LOG_OUTPUT_ON_FAILURE TRUE
                      BUILD_ALWAYS TRUE)
  add_library(haero STATIC IMPORTED GLOBAL)
  set_target_properties(haero PROPERTIES IMPORTED_LOCATION ${HAERO_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libhaero.a)

  # Build and install MAM4xx.
  set(MAM4XX_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/externals/mam4xx")
  set(MAM4XX_CMAKE_OPTS
    -DENABLE_TESTS=OFF
    -DENABLE_SKYWALKER=OFF
    -DCMAKE_INSTALL_PREFIX=${MAM4XX_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${mam_build_type}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DMAM4XX_HAERO_DIR=${HAERO_INSTALL_PREFIX}
    -DBUILD_SHARED_LIBS=OFF
    -DNUM_VERTICAL_LEVELS=${SCREAM_NUM_VERTICAL_LEV}
    )
  ExternalProject_Add(mam4xx_proj
                      PREFIX ${PROJECT_BINARY_DIR}/externals/mam4xx
                      SOURCE_DIR ${EXTERNALS_DIR}/mam4xx
                      BINARY_DIR ${PROJECT_BINARY_DIR}/externals/mam4xx
                      CMAKE_ARGS ${MAM4XX_CMAKE_OPTS}
                      DEPENDS haero_proj
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE
                      LOG_OUTPUT_ON_FAILURE TRUE
                      BUILD_ALWAYS TRUE)
  add_library(mam4xx STATIC IMPORTED GLOBAL)
  set_target_properties(mam4xx PROPERTIES IMPORTED_LOCATION ${MAM4XX_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libmam4xx.a)
  # Bring in MAM4xx-related targets by including the generated mam4xx.cmake
  #  list(APPEND CMAKE_MODULE_PATH ${MAM4XX_INSTALL_PREFIX}/share)
  #  include(mam4xx)

endif()
