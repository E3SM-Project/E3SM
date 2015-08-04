SET(SRC basic_tests.F90  driver.F90  global_vars.F90 ncdf_tests.F90 nc_set_log_level2.c)
ADD_DEFINITIONS(${PIO_DEFINITIONS})
FILE(COPY "./input.nl" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY "./not_netcdf.ieee" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Use to trace hdf5 errors with a special netcdf build
OPTION(HDF5_LOGGING "Turn on hdf5 logging (requires instrumented netcdf4)" OFF)
if(HDF5_LOGGING) 
  ADD_DEFINITIONS(-DLOGGING)
endif()

include(mpiexec)

INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS} )
LINK_DIRECTORIES(${PIO_LIB_DIR} )

# Default test timeout = 4 min (240s)
SET(DEFAULT_TEST_TIMEOUT 240)
ADD_EXECUTABLE(piotest ${SRC})
add_mpi_test("piotest" "${CMAKE_CURRENT_BINARY_DIR}/piotest" "" 4 ${DEFAULT_TEST_TIMEOUT})

if(${PIO_BUILD_TIMING}) 
  SET(TIMING_LINK_LIB timing)
endif()
TARGET_LINK_LIBRARIES(piotest piof pioc ${TIMING_LINK_LIB})

