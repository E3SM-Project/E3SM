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

# MAM4 aerosol support
if (SCREAM_ENABLE_MAM4)
  # We use CMake's ExternalProject capability to build and install Haero.
  include(ExternalProject)
  if (SCREAM_DOUBLE_PRECISION)
    set(HAERO_PRECISION "double")
  else()
    set(HAERO_PRECISION "single")
  endif()
  set(HAERO_CMAKE_OPTS
    -DHAERO_ENABLE_GPU=${Kokkos_ENAVBLE_CUDA}
    -DHAERO_ENABLE_MPI=ON
    -DHAERO_PRECISION=${HAERO_PRECISION}
    -DEKAT_INCLUDE_DIR=${PROJECT_BINARY_DIR}../../../externals/ekat/src
    -DEKAT_LIBRARY=${PROJECT_BINARY_DIR}/../../../externals/ekat/src/libekat.a)

  ExternalProject_Add(haero_proj
                      PREFIX ${PROJECT_BINARY_DIR}/externals/haero
                      SOURCE_DIR ../../../externals/haero
                      BINARY_DIR ${PROJECT_BINARY_DIR}/externals/haero
                      CMAKE_ARGS ${HAERO_CMAKE_OPTS}
                      DEPENDS ekat
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND ${MAKE} -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND ${MAKE} install
                      LOG_INSTALL TRUE
endif()
