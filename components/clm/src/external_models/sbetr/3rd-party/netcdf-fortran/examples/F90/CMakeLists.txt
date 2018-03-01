SET(example_F90_tests simple_xy_wr simple_xy_rd sfc_pres_temp_wr sfc_pres_temp_rd pres_temp_4D_wr pres_temp_4D_rd)
SET(example_F90_tests_source "")
FOREACH(F ${example_F90_tests}) 
  SET(example_F90_tests_source ${example_F90_test_source} ${F}.f90)
ENDFOREACH()

SET(prefix F90_tests)
FOREACH(F ${example_F90_tests})
  ADD_EXECUTABLE(${prefix}_${F} ${F}.f90)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff ${NETCDF_C_LIBRARY})
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F} 
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
	)
  ENDIF()
  ADD_TEST(${prefix}_${F} ${EXECUTABLE_OUTPUT_PATH}/${prefix}_${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${prefix}_${F} PROPERTY FOLDER "tests/")
  ENDIF()	
ENDFOREACH()

## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")

SET(CLEANFILES sfc_pres_temp.nc simple_xy.nc pres_temp_4D.nc simple_nc4.nc simple_xy_nc4.nc)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")
