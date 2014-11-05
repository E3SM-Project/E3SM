SET(SRC pioperformance.F90)
ADD_DEFINITIONS(${PIO_DEFINITIONS})

# Use to trace hdf5 errors with a special netcdf build
OPTION(HDF5_LOGGING "Turn on hdf5 logging (requires instrumented netcdf4)" OFF)
if(HDF5_LOGGING) 
  ADD_DEFINITIONS(-DLOGGING)
endif()

INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS} )
LINK_DIRECTORIES(${PIO_LIB_DIR} )



if(PIO_BUILD_TIMING MATCHES "ON") 
  ADD_EXECUTABLE(pioperf ${SRC})
  TARGET_LINK_LIBRARIES(pioperf piof pioc timing)
else()
  MESSAGE("Cannot build performance test without gptl timing library")
endif()

#ADD_TEST(pioperformance)
