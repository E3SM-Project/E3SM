# Omega Top-level CMakeLists.txt
#
# - supports two build modes:
#    1) standalone build
#    2) e3sm build


# Included scripts do automatic cmake_policy() PUSH and POP
cmake_policy(SET CMP0011 NEW)

# list command no longer ignores empty elements
cmake_policy(SET CMP0007 NEW)

# Only interpret if() arguments as variables or keywords when unquoted
cmake_policy(SET CMP0054 NEW)

#find_package() uses <PackageName>_ROOT variables.
cmake_policy(SET CMP0074 NEW)

###########################################################
# STEP 1: Setup                                           #
#                                                         #
# Assign values to the Omega build-controlling variables  #
###########################################################

# define variables and macros for Omega build
include(${CMAKE_CURRENT_SOURCE_DIR}/OmegaBuild.cmake)

# common pre-processing for standalone and e3sm mode
common()

# handle build modes
if (NOT DEFINED PROJECT_NAME)
  # enter standalone build

  cmake_minimum_required(VERSION 3.21) # the minimum version for HIP support

  # Collect machine and compiler info from CIME
  init_standalone_build()

  project(${OMEGA_PROJECT_NAME}
          LANGUAGES C CXX
         )

  set(CMAKE_CXX_STANDARD 17) # used in E3SM
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(LINKER_LANGUAGE C)     # needed to support Cray C compiler wrapper

  # update variables for standalone build
  setup_standalone_build()

else()
  # enter embedded build

  if(DEFINED E3SM_SOURCE_DIR) # E3SM_SOURCE_DIR = E3SM_ROOT/components

    # update variables for e3sm build
    setup_e3sm_build()

  else()
    message(FATAL_ERROR "UNKNOWN OMEGA BUILD MODE")

  endif()

endif()

###########################################################
# STEP 2: Update                                          #
#                                                         #
# update cmake & kokkos variables, and adjust               #
# build-controlling variables                             #
###########################################################

# update CMake and Kokkos variables from the build-controlling variables
update_variables()

###########################################
# place holder for further build-controls #
###########################################

# check if all-required variables have proper values
check_setup()


###########################################################
# STEP 3: Build                                           #
#                                                         #
# build Omega, and generate tests                          #
###########################################################

add_subdirectory(external)
add_subdirectory(src)

if(OMEGA_BUILD_TEST)

  enable_testing()
  add_subdirectory(test)

endif()

###########################################################
# STEP 4: Output                                          #
#                                                         #
# prepare outputs according to the output-controlling     #
# variables                                               #
###########################################################

wrap_outputs()
