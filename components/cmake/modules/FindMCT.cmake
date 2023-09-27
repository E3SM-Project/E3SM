# - Try to find MCT
#
# MCT does not provide a cmake config file, so we have to write our
# own find module.
#
# This can be controlled by setting the MCT_DIR (or, equivalently, the
# MCT_ROOT environment variable)
#
# Once done, this will define:
#
#   The "mct" target
#

# Look for libmct in INSTALL_SHAREDPATH/lib
find_library(MCT_LIB  mct  REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)
find_library(MPEU_LIB mpeu REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib $ENV{mct_ROOT})

# Create the interface library, and set target properties
add_library(mct INTERFACE)
target_link_libraries(mct INTERFACE ${MCT_LIB} ${MPEU_LIB})
target_include_directories(mct INTERFACE ${INSTALL_SHAREDPATH}/include)
target_compile_definitions(mct INTERFACE "MCT_INTERFACE")
