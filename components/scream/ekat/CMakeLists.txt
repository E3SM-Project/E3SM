cmake_minimum_required(VERSION 3.3)

# Add out own cmake goodies to cmake module search path
list(APPEND CMAKE_MODULE_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake
)

# EKAT has mostly C++, but some Fortran too
project (EKAT CXX Fortran)

###########################
###      EKAT TPLS      ###
###########################

option (EKAT_DISABLE_TPL_WARNINGS "Whether we should suppress warnings when compiling TPLs." OFF)

# WARNING: you CANNOT do list(APPEND var item1 ... item2) if var is a CACHE variable!
# Therefore, use an internal var during tpl parsing, then set a cache var ONCE at the end
set (EKAT_TPL_INCLUDE_DIRS_INTERNAL)
set (EKAT_TPL_LIBRARIES_INTERNAL)

# A good chunk of EKAT is a utility layer over kokkos
include(Kokkos)
list (APPEND EKAT_TPL_INCLUDE_DIRS_INTERNAL ${Kokkos_INCLUDE_DIRS})
list (APPEND EKAT_TPL_LIBRARIES_INTERNAL ${Kokkos_LIBRARIES})

# EKAT also has some yaml parsing utility
include(YamlCpp)
list (APPEND EKAT_TPL_INCLUDE_DIRS_INTERNAL ${YAMLCPP_INCLUDE_DIRS})
list (APPEND EKAT_TPL_LIBRARIES_INTERNAL ${YAMLCPP_LIBRARIES})

set (EKAT_TPL_LIBRARIES ${EKAT_TPL_LIBRARIES_INTERNAL} CACHE INTERNAL "List of EKAT's TPLs")
set (EKAT_TPL_INCLUDE_DIRS ${EKAT_TPL_INCLUDE_DIRS_INTERNAL} CACHE INTERNAL "List of EKAT's TPLs include directories")

############################
###  EKAT CONFIG OPTIONS ###
############################

# By default, we build ekat in double precision
option (EKAT_DOUBLE_PRECISION " Whether EKAT should be built with double precision floating points." ON)
option (EKAT_MPI_ERRORS_ARE_FATAL " Whether EKAT should crash when MPI errors happen." ON)
set (EKAT_PACK_SIZE       16 CACHE STRING
  " The number of scalars in a pack::Pack and Mask.
   Larger packs have good performance on conditional-free loops due to improved caching.")
set (EKAT_SMALL_PACK_SIZE  8 CACHE STRING
  " The number of scalars in a pack::SmallPack and SmallMask. Smaller packs can have better performance in loops with conditionals, since more of the packs will have masks with uniform value.")

# For some routines, SKX may have better performance with pksize=1
# The following lines are meant to set up another variable to store
# this pack size, which is usually equal to EKAT_PACK_SIZE,
# but may be equal to 1 if EKAT_POSSIBLY_NO_PACK is set to TRUE
# For example, by default we trigger this special behavior on SKX
if ("${KOKKOS_GMAKE_ARCH}" STREQUAL "SKX")
  set(DEFAULT_POSSIBLY_NO_PACK TRUE)
endif ()
set(EKAT_POSSIBLY_NO_PACK ${DEFAULT_POSSIBLY_NO_PACK} CACHE BOOL
  " Set possibly-no-pack to this value. You can set it to something else to restore packs on SKX for testing.")

set (DEFAULT_POSSIBLY_NO_PACK_SIZE ${EKAT_PACK_SIZE})
if (EKAT_POSSIBLY_NO_PACK)
  set (DEFAULT_POSSIBLY_NO_PACK_SIZE 1)
endif ()
set (EKAT_POSSIBLY_NO_PACK_SIZE ${DEFAULT_POSSIBLY_NO_PACK_SIZE})

############################
###     EKAT LIBRARY     ###
############################

# This is the only way to expose ekat's include dirs to host projects (if any)
set (EKAT_INCLUDE_DIRS_LIST
  ${EKAT_SOURCE_DIR}/src
  ${EKAT_BINARY_DIR}
  ${EKAT_BINARY_DIR}/src/ekat/ekat_f90_modules
  ${EKAT_TPL_INCLUDE_DIRS}
)

set (EKAT_INCLUDE_DIRS "${EKAT_INCLUDE_DIRS_LIST}" CACHE INTERNAL "EKAT library include dirs")

include(EkatSetCompilerFlags)
SetCompilerFlags()

add_subdirectory(src/ekat)

############################
###     EKAT TESTING     ###
############################

# By default, we DO build ekat tests. This may change in the future
option (EKAT_ENABLE_TESTS "Whether tests should be built." ON)

if (EKAT_ENABLE_TESTS)
  enable_testing()
  option (EKAT_STRICT_FP " Whether EKAT tests should adopt a strict fp model testing." ON)
  add_subdirectory(tests)
endif()

# Generate configuration file
include (EkatUtils)
EkatConfigFile (${EKAT_SOURCE_DIR}/ekat_config.h.in
                ${EKAT_BINARY_DIR}/ekat_config.h
                ${EKAT_BINARY_DIR}/ekat_config.f)
