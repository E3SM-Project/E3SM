#==============================================================================
#
#  API documentation with Doxygen
#
#==============================================================================

find_package(Doxygen)

# This supports the build with/witout async code. Once async code is
# fully merged, remove the definition of C_SRC_FILES and its mention
# in Doxyfile.in for simplicity.
if (PIO_ENABLE_ASYNC)
  SET(C_SRC_FILES "@CMAKE_CURRENT_SOURCE_DIR@/../src/clib/bget.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pioc.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pioc_sc.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_darray_async.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_get_nc_async.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_internal.h @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_nc4.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_put_nc_async.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_spmd.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/bget.h @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pioc_support.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_lists.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_nc_async.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_varm.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/dtypes.h  @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_file.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio.h @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_msg.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_rearrange.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/topology.c")
else ()
  SET(C_SRC_FILES "@CMAKE_CURRENT_SOURCE_DIR@/../src/clib/bget.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pioc.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pioc_sc.c  @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_internal.h @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_nc4.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_spmd.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/bget.h @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pioc_support.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_darray.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_get_nc.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_lists.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_put_nc.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_varm.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/dtypes.h @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_file.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio.h @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_msg.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_nc.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/pio_rearrange.c @CMAKE_CURRENT_SOURCE_DIR@/../src/clib/topology.c")
endif ()

if(DOXYGEN_FOUND)
  # Process the Doxyfile using options set during configure.
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

  # Copy necessary files.
  add_custom_target(doc
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/customdoxygen.css 
    ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml 
    ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.sty 
    ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)
