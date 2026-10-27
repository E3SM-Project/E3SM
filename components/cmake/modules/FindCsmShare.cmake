# - Try to find CsmShare
#
# CsmShare does not provide a cmake config file, so we have to write our
# own find module.
#
# This can be controlled by setting the CsmShare_DIR (or, equivalently, the
# csm_share_ROOT environment variable)
#
# Once done, this will define:
#
#   The "csm_share" target
#

if (TARGET csm_share)
  return()
endif()

# Build the name of the path where libcsm_share should be located
if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()
set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")

# Look for libcsm_share in the complex path we built above
find_library(CSM_SHARE_LIB csm_share REQUIRED HINTS ${CSM_SHARE})

set(CSM_LIBS  ${CSM_SHARE_LIB};mct)
set(CSM_INCLS ${CSM_SHARE})
if (USE_MOAB)
  list(APPEND CSM_LIBS ${MOAB_LIBRARIES})
  list(APPEND CSM_INCLS ${MOAB_INCLUDE_DIRS})
endif()

# Grab mct and spio if we don't already have it
find_package(MCT REQUIRED)
if (NOT DEFINED ENV{SCORPIO_ROOT})
  set(ENV{SCORPIO_ROOT} ${INSTALL_SHAREDPATH})
endif()

#set(CMAKE_FIND_DEBUG_MODE TRUE)
if (NOT SCORPIO_FOUND)
  find_package(SCORPIO REQUIRED)
endif()

# Create the interface library, and set target properties
add_library (csm_share INTERFACE)
target_link_libraries (csm_share INTERFACE ${CSM_LIBS})
target_include_directories(csm_share INTERFACE ${CSM_INCLS})

target_link_libraries(csm_share INTERFACE scorpio::piof scorpio::pioc)
