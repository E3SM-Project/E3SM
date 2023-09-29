# - Try to find PIO
#
# Placeholder code until spio gets a proper config.cmake
#
# Once done, this will define:
#
#   The "spio" target
#

if (TARGET spio)
  return()
endif()

if (NOT PIO_LIBDIR)
  set(PIO_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (PIO_VERSION STREQUAL 2)
  # This is a pio2 library
  set(PIOLIBS "${PIO_LIBDIR}/libpiof.a;${PIO_LIBDIR}/libpioc.a")
else()
  # This is a pio1 library
  set(PIOLIBS "${PIO_LIBDIR}/libpio.a")
endif()

find_package(NETCDF REQUIRED)

# Create the interface library, and set target properties
add_library(spio INTERFACE)
target_link_libraries(spio INTERFACE ${PIOLIBS};netcdf)
