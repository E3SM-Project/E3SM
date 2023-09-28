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

# Build the name of the path where libcsm_share should be located
if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()
set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")

# Look for libcsm_share in the complex path we built above
find_library(CSM_SHARE_LIB csm_share REQUIRED PATHS ${CSM_SHARE})

# Grab mct if we don't already have it
if (NOT TARGET mct)
  find_package(MCT REQUIRED)
endif()

# Create the interface library, and set target properties
add_library (csm_share INTERFACE)
target_link_libraries (csm_share INTERFACE ${CSM_SHARE_LIB};mct;${PIOLIBS})
target_include_directories(csm_share INTERFACE ${CSM_SHARE})

# Link against piof. Don't worry about this for now. Fix once spio
# has a proper config.cmake.
# target_link_libraries(csm_share INTERFACE piof)
