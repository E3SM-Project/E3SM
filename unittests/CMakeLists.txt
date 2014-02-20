SET(SRC basic_tests.F90  driver.F90  global_vars.F90 ncdf_tests.F90 nc_set_log_level2.c)
ADD_DEFINITIONS(${PIO_DEFINITIONS})

# Use to trace hdf5 errors with a special netcdf build
OPTION(HDF5_LOGGING "Turn on hdf5 logging (requires instrumented netcdf4)" OFF)
if(HDF5_LOGGING) 
  ADD_DEFINITIONS(-DLOGGING)
endif()

INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS} )
LINK_DIRECTORIES(${PIO_LIB_DIR} )

ADD_EXECUTABLE(piotest ${SRC})

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../timing)
  TARGET_LINK_LIBRARIES(piotest pio timing)
ELSE()
  TARGET_LINK_LIBRARIES(piotest pio)
ENDIF()
ADD_TEST(piotest piotest)
