cmake_minimum_required(VERSION 2.8)
project(clm45_tests Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CESM_CMAKE_MODULE_DIRECTORY})
include(CESM_utils)

set(CLM_ROOT "..")
set(CESM_ROOT "${CLM_ROOT}/../../..")

# Add source directories from other models (csm_share, etc.). This should be
# done first, so that in case of name collisions, the CLM versions take
# precedence (when there are two files with the same name, the one added later
# wins).
add_subdirectory(${CESM_ROOT}/models/csm_share/shr csm_share)

# Add CLM source directories (these add their own test directories)
add_subdirectory(${CLM_ROOT}/src/utils clm_utils)
add_subdirectory(${CLM_ROOT}/src/biogeophys clm_biogeophys)
add_subdirectory(${CLM_ROOT}/src/dyn_subgrid clm_dyn_subgrid)
add_subdirectory(${CLM_ROOT}/src/main clm_main)

# Add general unit test directories (mocked out files, etc.)
add_subdirectory(unit_test_mocks)
add_subdirectory(unit_test_shr)

# Remove shr_sys_mod from share_sources.
# This is needed because we want to use the mock shr_sys_mod in place of the real one
#
# TODO: this should be moved into a general-purpose function in Sourcelist_utils.
# Then this block of code could be replaced with a single call, like:
# remove_source_file(${share_sources} "shr_sys_mod.F90")
foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_sys_mod.F90" match_found ${sourcefile})
  if(match_found)
    list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()

# TODO: Figure out where to put declare_generated_dependencies, or whether to just leave that out

# Build libraries containing stuff needed for the unit tests. 
# Eventually, these add_library calls should probably be distributed into the correct location, rather than being in this top-level CMakeLists.txt file.
add_library(csm_share ${share_sources})
add_library(clm ${clm_sources})
add_dependencies(clm csm_share)

# We need to look for header files here, in order to pick up shr_assert.h
include_directories(${CESM_ROOT}/models/csm_share/shr)

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add the test directories
# Note: it's possible that these could be added by each source directory that
# has tests in it. However, it appears that the order needs to be done
# carefully: for example, include_directories and link_directories needs to be
# done before adding the tests themselves.
add_subdirectory(${CLM_ROOT}/src/biogeophys/test clm_biogeophys_test)
add_subdirectory(${CLM_ROOT}/src/dyn_subgrid/test clm_dyn_subgrid_test)
add_subdirectory(${CLM_ROOT}/src/main/test clm_main_test)
