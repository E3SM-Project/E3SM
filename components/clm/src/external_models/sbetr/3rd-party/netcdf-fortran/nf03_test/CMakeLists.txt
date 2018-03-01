SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_BUILD_TYPE "RelWithDebInfo")
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

SET(TESTS)

# Need some stuff from ../nf_test
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/nf_test ${NC4F_BINARY_DIR}/libsrc
   ${NC4F_BINARY_DIR}/fortran)
# All tests need to link to fortran and C libraries.
LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/fortran)
#FIND_LIBRARY(NETCDF_C_LIBRARY NAMES netcdf libnetcdf)
#IF(NOT NETCDF_C_LIBRARY)
#  MESSAGE(FATAL_ERROR "Cannot find netCDF C library, needed for netCDF Fortran tests.")
#ENDIF()
#FIND_PATH(NETCDF_C_INCLUDE_DIR netcdf.h)

# nf03_test is the main test program.
SET(check_PROGRAMS nf03_test)
LIST(APPEND TESTS nf03_test)
SET(nf03_test_SOURCES module_tests.F90 f03lib_f_interfaces.f90 test03_get.F test03_put.F nf03_error.F
  nf03_test.F test03_read.F test03_write.F util03.F f03lib.c)

IF(USE_NETCDF4)
  LIST(APPEND check_PROGRAMS f03tst_vars f03tst_vars2 f03tst_vars3 f03tst_vars4 f03tst_vars5
    f03tst_vars6 f03tst_types f03tst_types2 f03tst_types3 f03tst_groups)
  SET(f03tst_types_SOURCES module_tests.F90 f03tst_types.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_types2_SOURCES module_tests.F90 f03tst_types2.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_types3_SOURCES module_tests.F90 f03tst_types3.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_vars_SOURCES module_tests.F90 f03tst_vars.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_vars2_SOURCES module_tests.F90 f03tst_vars2.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_vars3_SOURCES module_tests.F90 f03tst_vars3.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_vars4_SOURCES module_tests.F90 f03tst_vars4.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_vars5_SOURCES module_tests.F90 f03tst_vars5.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_vars6_SOURCES module_tests.F90 f03tst_vars6.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  SET(f03tst_groups_SOURCES module_tests.F90 f03tst_groups.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  LIST(APPEND TESTS f03tst_vars f03tst_vars2 f03tst_vars3 f03tst_vars4 f03tst_vars5
    f03tst_vars6 f03tst_types f03tst_types2 f03tst_types3 f03tst_groups)
  LIST(APPEND CLEANFILES ftst_vars.nc ftst_vars2.nc ftst_vars3.nc ftst_vars4.nc
    ftst_vars5.nc ftst_vars6.nc ftst_types.nc ftst_types2.nc ftst_types3.nc ftst_groups.nc)

  # This is a netCDF-4 V2 test program.
  IF(BUILD_V2)
    LIST(APPEND check_PROGRAMS f03tst_v2)
    LIST(APPEND TESTS f03tst_v2)
    LIST(APPEND CLEANFILES ftst_v2.nc)
    SET(f03tst_v2_SOURCES f03tst_v2.F)
  ENDIF(BUILD_V2)

  # Test parallel I/O
  IF(TEST_PARALLEL)
    LIST(APPEND check_PROGRAMS f03tst_parallel)
    LIST(APPEND CLEANFILES ftst_parallel.nc)
    SET(f03tst_parallel_SOURCES f03tst_parallel.F)
  ENDIF(TEST_PARALLEL)
ENDIF(USE_NETCDF4)

# This is the fortran v2 test. It depends on the utilities being built
# to generate it's input file.
IF(BUILD_V2)
  LIST(APPEND check_PROGRAMS f03test)
  LIST(APPEND TESTS f03test)
  LIST(APPEND CLEANFILES test.nc)
  SET(f03test_SOURCES module_tests.F90 f03test.F f03lib.c)
ENDIF(BUILD_V2)

LIST(APPEND CLEANFILES tests.mod)

# test03_get.F and test03_put.f need to be distributed, so that the user
# need not have m4.
SET(EXTRA_DIST test03_get.F test03_put.F fills.cdl create_fills.sh run_f77_par_test.sh)

# Did the user build the V2 F77 API? If so, run this test.
IF(BUILD_V2)
  LIST(APPEND check_PROGRAMS tst03_f77_v2)
  LIST(APPEND TESTS tst03_f77_v2)
  LIST(APPEND CLEANFILES tst_f77_v2.nc)
  SET(tst03_f77_v2_SOURCES module_tests.F90 tst03_f77_v2.F)
ENDIF(BUILD_V2)

# Need a copy of ref_fills.nc for ftest
execute_process(COMMAND cp ${CMAKE_SOURCE_DIR}/nf_test/ref_fills.nc ${CMAKE_CURRENT_BINARY_DIR}/fills.nc)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/fortran ${CMAKE_SOURCE_DIR}/libsrc ${CMAKE_SOURCE_DIR}
  ${CMAKE_CPP_FLAGS} ${NETCDF_C_INCLUDE_DIR})
LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/fortran ${CMAKE_SOURCE_DIR}/libsrc)
FOREACH(t ${TESTS})
  ADD_EXECUTABLE(${t} ${${t}_SOURCES})
  TARGET_LINK_LIBRARIES(${t} netcdff)
  ADD_TEST(${t} ${t})
ENDFOREACH()

# FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.sh ${CMAKE_CURRENT_SOURCE_DIR}/*.cdl)
# FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE)

# FOREACH(F ${CDL_EXAMPLE_TESTS})
# 	add_sh_test(cdl ${F})
# ENDFOREACH()

# SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")
