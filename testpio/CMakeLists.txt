SET(TEMPSRCF90 alloc_mod.F90)
FOREACH(tempfile IN LISTS TEMPSRCF90)
ADD_CUSTOM_COMMAND(
	OUTPUT ${tempfile}
	COMMAND ${GENF90_PATH}/genf90.pl ${CMAKE_CURRENT_SOURCE_DIR}/${tempfile}.in > ${tempfile}
	DEPENDS  ${tempfile}.in
)
ENDFOREACH()

SET(SRC check_mod.F90  gdecomp_mod.F90  kinds_mod.F90  namelist_mod.F90  
            testpio.F90  utils_mod.F90 ${TEMPSRCF90})
SET(WSSRC wstest.c)

INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})
LINK_DIRECTORIES(${PIO_LIB_DIR})
ADD_EXECUTABLE(testpio ${SRC})
ADD_EXECUTABLE(wstest ${WSSRC})
if(PIO_BUILD_TIMING MATCHES "ON") 
  TARGET_LINK_LIBRARIES(testpio timing  )
  TARGET_LINK_LIBRARIES(wstest timing  )
endif()
if(NETCDF_FOUND)
  TARGET_LINK_LIBRARIES(testpio piof pioc ${NETCDF_Fortran_LIBRARY})
else()
  TARGET_LINK_LIBRARIES(testpio piof pioc)
endif()
TARGET_LINK_LIBRARIES(wstest pioc)





