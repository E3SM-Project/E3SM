cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

project(drv_tests Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CESM_CMAKE_MODULE_DIRECTORY})
include(CESM_utils)

set(DRV_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(CESM_ROOT "${DRV_ROOT}/../..")

# ------------------------------------------------------------------------
# Build mct
# TODO: move this code to build mct into a standalone utility somewhere
# ------------------------------------------------------------------------
set(MCT_ROOT "${CESM_ROOT}/utils/mct")
ExternalProject_add(mct_project
  SOURCE_DIR ${MCT_ROOT}
  CONFIGURE_COMMAND ${MCT_ROOT}/configure --enable-mpiserial --enable-debugging --prefix=${CMAKE_CURRENT_BINARY_DIR}
  BUILD_COMMAND make
  )
# This copy_makefiles step is needed because mct currently doesn't support an
# out-of-source build. I am replicating what is done for the CESM system build.
ExternalProject_add_step(mct_project copy_makefiles
  DEPENDEES configure
  DEPENDERS build
  WORKING_DIRECTORY <BINARY_DIR>
  COMMAND cp -p <SOURCE_DIR>/Makefile .
  COMMAND mkdir -p mct
  COMMAND cp -p <SOURCE_DIR>/mct/Makefile mct/
  COMMAND mkdir -p mpeu
  COMMAND cp -p <SOURCE_DIR>/mpeu/Makefile mpeu/
  COMMAND mkdir -p mpi-serial
  COMMAND cp -p <SOURCE_DIR>/mpi-serial/Makefile mpi-serial/
  )
# This copy_files_for_install step seems to be needed in order for MCT's 'make
# install' to work properly with this out-of-source build. It is not needed by
# the CESM system build, because (it seems) the system build replaces 'make
# install' with its own copies.
ExternalProject_add_step(mct_project copy_files_for_install
  DEPENDEES copy_makefiles
  DEPENDERS install
  WORKING_DIRECTORY <BINARY_DIR>
  COMMAND cp -p <SOURCE_DIR>/mkinstalldirs .
  COMMAND cp -p <SOURCE_DIR>/install-sh .
  # Copying *.h doesn't work, so I'm copying each needed header file
  # individually
  COMMAND cp -p <SOURCE_DIR>/mpi-serial/mpif.real4double8.h mpi-serial/
  COMMAND cp -p <SOURCE_DIR>/mpi-serial/mpi.h mpi-serial/
  )
# ------------------------------------------------------------------------
# Done MCT build
# ------------------------------------------------------------------------

add_definitions(
  -DNUM_COMP_INST_ATM=1
  -DNUM_COMP_INST_LND=1
  -DNUM_COMP_INST_OCN=1
  -DNUM_COMP_INST_ICE=1
  -DNUM_COMP_INST_GLC=1
  -DNUM_COMP_INST_WAV=1
  -DNUM_COMP_INST_ROF=1
  )

# The following definitions are needed when building with the mpi-serial library
add_definitions(-DNO_MPI2 -DNO_MPIMOD)

# Add source directories from other models (csm_share, etc.). This should be
# done first, so that in the case of name collisions, the drv versions take
# precedence (when there are two files with the same name, teh one added later
# wins).
add_subdirectory(${CESM_ROOT}/cime/share/csm_share/shr csm_share)
add_subdirectory(${CESM_ROOT}/cime/share/csm_share/unit_test_stubs/pio pio)
add_subdirectory(${CESM_ROOT}/cime/share//esmf_wrf_timemgr esmf_wrf_timemgr)

# Add drv source directories
add_subdirectory(${DRV_ROOT}/shr drv_shr)
add_subdirectory(${DRV_ROOT}/driver drv_driver)

# Add general unit test directories (stubbed out files, etc.)
add_subdirectory(utils drv_unit_test_utils)

# TODO: Figure out where to put declare_generated_dependencies, or whether to
# just leave that out

# Build libraries containing stuff needed for the unit tests. 
# Eventually, these add_library calls should probably be distributed into the
# correct location, rather than being in this top-level CMakeLists.txt file.
# Note that we are including the stub pio in the csm_share library for simplicity.
add_library(csm_share ${share_sources};${share_mct_sources};${share_pio_sources};${pio_sources})
add_dependencies(csm_share mct_project)
add_library(esmf_wrf_timemgr ${esmf_wrf_timemgr_sources})
add_dependencies(esmf_wrf_timemgr csm_share)
add_library(drv ${drv_sources})
add_dependencies(drv csm_share esmf_wrf_timemgr)

# We need to look for header files here, in order to pick up shr_assert.h
include_directories(${CESM_ROOT}/cime/share/csm_share/shr)

# And we need to look for header files here, for some include files needed by
# the esmf_wrf_timemgr code
include_directories(${CESM_ROOT}/cime/share/esmf_wrf_timemgr)

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
link_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

# Add the test directories
add_subdirectory(avect_wrapper_test)
add_subdirectory(seq_map_test)
