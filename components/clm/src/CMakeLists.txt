# AGS: CMake logic already existed -- preserving that functionalituy
#       within !ACME_BUILD.
if (NOT ACME_BUILD)

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
add_subdirectory(${CLM_ROOT}/src/external_models/fates/main ed_main)
add_subdirectory(${CLM_ROOT}/src/external_models/fates/biogeophys ed_biogeophys)
add_subdirectory(${CLM_ROOT}/src/external_models/fates/biogeochem ed_biogeochem)
add_subdirectory(${CLM_ROOT}/src/external_models/fates/fire ed_fire)


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

ELSE() #ACME_BUILD
  MESSAGE("   + CLM:")


  IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/main/ncdio_pio.F90)
  #IF (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/shr_infnan_mod.F90)
    # Do genf90 on 2 input files
    SET (GENF90 ${ACME_SOURCE_DIR}/cime/src/externals/genf90/genf90.pl)
    MESSAGE("     - Generating ncdio_pio.mod ${GENF90}")
    execute_process (COMMAND ${GENF90}
                     ${CMAKE_CURRENT_SOURCE_DIR}/main/ncdio_pio.F90.in
                     OUTPUT_FILE ncdio_pio.F90
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/main)
  ENDIF()
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "${CMAKE_CURRENT_SOURCE_DIR}/main/ncdio_pio.F90;${CMAKE_CURRENT_SOURCE_DIR}/main/dtypes.h")


  #AGS: Globbing all clm files in one library. There appear to
  #     be to many inter-directory dependencies to tease out
  #     separate libraries with one-way dependencies.
  FILE(GLOB CLM_MAIN_SOURCE         main/*F90 ED/main/*F90)
  FILE(GLOB CLM_UTILS_SOURCE        utils/*F90)
  FILE(GLOB CLM_EDFIRE_SOURCE       ED/fire/*F90)
  FILE(GLOB CLM_BIOGEOPHYS_SOURCE   biogeophys/*F90 ED/biogeophys/*F90)
  FILE(GLOB CLM_BIOGEOCHEM_SOURCE   biogeochem/*F90 ED/biogeochem/*F90)
  FILE(GLOB CLM_DYN_SUBGRID_SOURCE  dyn_subgrid/*F90)
  FILE(GLOB CLM_BETR_SOURCE         betr/*F90 betr/betr_core/*F90
         betr/betr_math/*F90 betr/bgc_century/*F90 betr/bgc_sminn/*F90)
  FILE(GLOB CLM_CPL_SOURCE  cpl/*F90)

# Look into CMake ExternalProject for FATES
  FILE(GLOB CLM_FATES_SOURCE 
         external_models/fates/main/*F90
         external_models/fates/biogeophys/*F90
         external_models/fates/biogeochem/*F90
         external_models/fates/fire/*F90
         external_models/emi/src/*F90
         external_models/mpp/src/mpp/util/*F90)

  message("   clm build: add_definitions: NDEBUG Required!!")
#  add_definitions( -mcmodel=medium -fconvert=big-endian
#-ffree-line-length-none -ffixed-line-length-none  -O2 -DLINUX  -DNDEBUG
#-DMCT_INTERFACE -DNO_R16 -DHAVE_SLASHPROC -ffree-form)
  add_definitions(-O2 -DNDEBUG)

  INCLUDE_DIRECTORIES (
          ${SHR_BINARY_DIR}
          ${GPTL_BINARY_DIR}  #perf_mod
          ${ESMF_BINARY_DIR} # esmf.mod
          ${ACME_SHR_DIR}  #shr_assert.h
          ${SEQ_BINARY_DIR}
  )

  #  ${MPI_Fortran_INCLUDE_PATH}
  ADD_LIBRARY (clm
               ${CLM_MAIN_SOURCE} ${CLM_UTILS_SOURCE}
               ${CLM_BIOGEOPHYS_SOURCE} ${CLM_BIOGEOCHEM_SOURCE}
               ${CLM_EDFIRE_SOURCE} ${CLM_DYN_SUBGRID_SOURCE}
               ${CLM_BETR_SOURCE} ${CLM_CPL_SOURCE}
               ${CLM_FATES_SOURCE}
              )
  ADD_DEPENDENCIES(clm shrsys seq shr gptl)

ENDIF()
