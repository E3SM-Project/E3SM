# Omega Top-level CMakeLists.txt
#
# - supports two build modes:
#    1) standalone build
#    2) e3sm build


###########################################################
# STEP 1: Setup                                           #
#                                                         #
# Assign values to the Omega build-controlling variables  #
###########################################################

# define variables and macros for Omega build
include(${CMAKE_CURRENT_SOURCE_DIR}/OmegaBuild.cmake)


# handle build modes 
if (NOT DEFINED PROJECT_NAME)
  # enter standalone build

  preset()

  cmake_minimum_required(VERSION 3.16) # used in E3SM

  project(${OMEGA_PROJECT_NAME}
             LANGUAGES CXX
         )

  set(CMAKE_CXX_STANDARD 17) # used in E3SM

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
# update cmake & yakl variables, and adjust               #
# build-controlling variables                             #
###########################################################

# update CMake and YAKL variables from the build-controlling variables
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
