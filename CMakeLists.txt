cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

project(cime_tests Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CESM_CMAKE_MODULE_DIRECTORY})
include(CESM_utils)
find_package(NETCDF)
include_directories(${NETCDF_INCLUDE_DIRS})

set(CIME_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# TODO: Some of the below should be done in the relevant directories, not in
# this top level CMakeLists.

# ------------------------------------------------------------------------
# Build mct
# ------------------------------------------------------------------------
set(MCT_ROOT "${CIME_ROOT}/externals/mct")

ExternalProject_add(mct_project
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}
  SOURCE_DIR ${MCT_ROOT}
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/mct
  CONFIGURE_COMMAND ${MCT_ROOT}/configure --enable-debugging --prefix=${CMAKE_CURRENT_BINARY_DIR}
  BUILD_COMMAND make SRCDIR=${MCT_ROOT}
  # Leave things in <BINARY_DIR> rather than "installing", because we have
  # no need to move things around inside of the CMake binary directory.
  INSTALL_COMMAND :
  )
# This copy_makefiles step is needed because mct currently doesn't support an
# out-of-source build. I am replicating what is done for the CIME system build.
ExternalProject_add_step(mct_project copy_makefiles
  DEPENDEES configure
  DEPENDERS build
  WORKING_DIRECTORY <BINARY_DIR>
  COMMAND cp -p <SOURCE_DIR>/Makefile .
  COMMAND mkdir -p mct
  COMMAND cp -p <SOURCE_DIR>/mct/Makefile mct/
  COMMAND mkdir -p mpeu
  COMMAND cp -p <SOURCE_DIR>/mpeu/Makefile mpeu/
  )

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR}/mct/mct)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/mct/mpeu)
link_directories(${CMAKE_CURRENT_BINARY_DIR}/mct/mct)
link_directories(${CMAKE_CURRENT_BINARY_DIR}/mct/mpeu)

# ------------------------------------------------------------------------
# Done MCT build
# ------------------------------------------------------------------------

# Now a bunch of includes for share code.

# csm_share (we don't build it here because it seems to be built differently
# by different tests?)
set(SHARE_ROOT "${CIME_ROOT}/share/csm_share")
add_subdirectory(${SHARE_ROOT}/shr csm_share)
include_directories(${SHARE_ROOT}/include)

# esmf_wrf_timemgr not built here because it depends on csm_share.
add_subdirectory(${CIME_ROOT}/share/esmf_wrf_timemgr esmf_wrf_timemgr)
include_directories(${CIME_ROOT}/share/esmf_wrf_timemgr)

# Now the actual test directories.
add_subdirectory(driver_cpl/unit_test)
add_subdirectory(share/csm_share/test/unit)
