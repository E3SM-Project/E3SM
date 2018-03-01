CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_BUILD_TYPE "RelWithDebInfo")
SET(CMAKE_INCLUDE_CURRENT_DIR ON)


# Process these files with m4.
SET(m4_SOURCES test_get test_put)
foreach (f ${m4_SOURCES})
		   GEN_m4(${f})
endforeach(f)


IF (BUILD_F03)
  SET(nf_test_SOURCES
    f03lib_f_interfaces.f90 test_get.F test_put.F nf_error.F nf_test.F
    test_read.F test_write.F util.F f03lib.c tests.inc)
ELSE()
  SET(nf_test_SOURCES
    test_get.F test_put.F nf_error.F nf_test.F
    test_read.F test_write.F util.F fortlib.c tests.inc
    )
ENDIF()

SET(TESTS nf_test)
IF (BUILD_V2)
  SET(check_PROGRAMS tst_f77_v2)
  SET(tst_f77_v2_SOURCES tst_f77_v2.F)
  SET(TESTS ${TESTS} tst_f77_v2)
ENDIF(BUILD_V2)

# Is the user building netCDF-4?
if (USE_NETCDF4)
  # Add these netCDF-4 test programs.
  SET(nc4_check_PROGRAMS
    ftst_vars ftst_vars2 ftst_vars3 ftst_vars4 ftst_vars5
    ftst_vars6 ftst_types ftst_types2 ftst_types3 ftst_groups ftst_path ftst_rengrps)
  SET(check_PROGRAMS ${check_PROGRAMS} ${nc4_check_PROGRAMS})
  SET(TESTS ${TESTS} ${nc4_check_PROGRAMS})

  if (BUILD_F03)
    SET(ftst_types_SOURCES ftst_types.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_types2_SOURCES ftst_types2.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_types3_SOURCES ftst_types3.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_vars_SOURCES ftst_vars.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_vars2_SOURCES ftst_vars2.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_vars3_SOURCES ftst_vars3.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_vars4_SOURCES ftst_vars4.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_vars5_SOURCES ftst_vars5.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_vars6_SOURCES ftst_vars6.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_groups_SOURCES ftst_groups.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_path_SOURCES ftst_path.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
    SET(ftst_rengrps_SOURCES ftst_rengrps.F f03lib_f_interfaces.f90 f03lib.c handle_err.F)
  else()
    SET(ftst_types_SOURCES ftst_types.F fortlib.c handle_err.F)
    SET(ftst_types2_SOURCES ftst_types2.F fortlib.c handle_err.F)
    SET(ftst_types3_SOURCES ftst_types3.F fortlib.c handle_err.F)
    SET(ftst_vars_SOURCES ftst_vars.F fortlib.c handle_err.F)
    SET(ftst_vars2_SOURCES ftst_vars2.F fortlib.c handle_err.F)
    SET(ftst_vars3_SOURCES ftst_vars3.F fortlib.c handle_err.F)
    SET(ftst_vars4_SOURCES ftst_vars4.F fortlib.c handle_err.F)
    SET(ftst_vars5_SOURCES ftst_vars5.F fortlib.c handle_err.F)
    SET(ftst_vars6_SOURCES ftst_vars6.F fortlib.c handle_err.F)
    SET(ftst_groups_SOURCES ftst_groups.F fortlib.c handle_err.F)
    SET(ftst_path_SOURCES ftst_path.F fortlib.c handle_err.F)
    SET(ftst_rengrps_SOURCES ftst_rengrps.F fortlib.c handle_err.F)
  endif(BUILD_F03)
  # Add these netCDF-4 f90 test programs.
  SET(nc4_f90_PROGRAMS
    f90tst_vars tst_types tst_types2 f90tst_vars_vlen
    tst_f90_nc4 f90tst_grps f90tst_fill f90tst_fill2 f90tst_vars3
    f90tst_vars4 f90tst_vars2)
  SET(check_PROGRAMS ${check_PROGRAMS} ${nc4_f90_PROGRAMS})
  SET(TESTS ${TESTS} ${nc4_f90_PROGRAMS})
  SET(f90tst_vars_SOURCES f90tst_vars.f90)
  SET(tst_types_SOURCES tst_types.f90)
  SET(tst_types2_SOURCES tst_types2.f90)
  SET(tst_f90_nc4_SOURCES tst_f90_nc4.f90)
  SET(f90tst_vars_vlen_SOURCES f90tst_vars_vlen.f90)
  SET(f90tst_grps_SOURCES f90tst_grps.f90)
  SET(f90tst_fill_SOURCES f90tst_fill.f90)
  SET(f90tst_fill2_SOURCES f90tst_fill2.f90)
  SET(f90tst_vars2_SOURCES f90tst_vars2.f90)
  SET(f90tst_vars3_SOURCES f90tst_vars3.f90)
  SET(f90tst_vars4_SOURCES f90tst_vars4.f90)

  # This is a netCDF-4 V2 test program.
  IF (BUILD_V2)
    SET(nc4_v2_PROGRAMS ftst_v2)
    SET(check_PROGRAMS ${check_PROGRAMS} ${nc4_v2_PROGRAMS})
    SET(TESTS ${TESTS} ${nc4_v2_PROGRAMS})
    SET(ftst_v2_SOURCES ftst_v2.F)
  ENDIF (BUILD_V2)

  # This is the netCDF-4 F90 large file test.
  IF (LARGE_FILE_TESTS)
    SET(nc4_largefile_PROGRAMS tst_flarge)
    SET(check_PROGRAMS ${check_PROGRAMS} ${nc4_largefile_PROGRAMS})
    SET(TESTS ${TESTS} ${nc4_largefile_PROGRAMS})
    SET(tst_flarge_SOURCES tst_flarge.f90)
  ENDIF(LARGE_FILE_TESTS)

  # This is an f90 benchmark.
  IF(BUILD_BENCHMARKS)
    SET(nc4_io_PROGRAMS tst_io)
    SET(check_PROGRAMS ${check_PROGRAMS} ${nc4_io_PROGRAMS})
    SET(TESTS ${TESTS} ${nc4_io_PROGRAMS})
    SET(tst_io_SOURCES tst_io.f90)
  ENDIF(BUILD_BENCHMARKS)

  # Test parallel I/O.
  IF (TEST_PARALLEL)
    SET(check_PROGRAMS ${check_PROGRAMS} f90tst_parallel f90tst_parallel2 f90tst_parallel3
      f90tst_nc4_par)

    build_bin_test(f90tst_parallel ".f90")
    build_bin_test(f90tst_parallel2 ".f90")
    build_bin_test(f90tst_parallel3 ".f90")
    build_bin_test(f90tst_nc4_par ".f90")

    #SET(f90tst_parallel_fill_SOURCES f90tst_parallel_fill.f90)

  ENDIF (TEST_PARALLEL)

  # Test parallel I/O for F77.
  IF (TEST_PARALLEL)
    SET(check_PROGRAMS ${check_PROGRAMS} ftst_parallel ftst_parallel_nasa)
    SET(SCRIPT_TESTS ${SCRIPT_TESTS} run_f77_par_test)

    # Add executables for f77_par_tst
    build_bin_test(ftst_parallel ".F")
    build_bin_test(ftst_parallel_nasa ".F")

  ENDIF (TEST_PARALLEL)

endif(USE_NETCDF4)

# This is the fortran v2 test. It depends on the utilities being built
# to generate its input file.
IF (BUILD_V2)
  SET(check_PROGRAMS ${check_PROGRAMS} ftest)
  SET(TESTS ${TESTS} ftest)
  IF (BUILD_F03)
    SET(ftest_SOURCES ftest.F f03lib.c)
  ELSE()
    SET(ftest_SOURCES ftest.F fortlib.c)
  ENDIF(BUILD_F03)
ENDIF(BUILD_V2)

SET(check_PROGRAMS ${check_PROGRAMS} tst_f90)
SET(TESTS ${TESTS} tst_f90)
SET(tst_f90_SOURCES tst_f90.f90)

FIND_LIBRARY(NETCDF_C_LIBRARY NAMES netcdf libnetcdf)
FIND_PATH(NETCDF_C_INCLUDE_DIR netcdf.h)
IF(NOT NETCDF_C_LIBRARY)
  MESSAGE(FATAL_ERROR "Cannot find netCDF C library, needed for netCDF Fortran tests.")
ENDIF()

# Need a copy of ref_fills.nc for ftest
execute_process(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/ref_fills.nc
  ${CMAKE_CURRENT_BINARY_DIR}/fills.nc)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/fortran ${CMAKE_SOURCE_DIR}/libsrc
  ${CMAKE_SOURCE_DIR} ${CMAKE_CPP_FLAGS}
  ${NETCDF_C_INCLUDE_DIR} ${NC4F_BINARY_DIR}/fortran
  ${NC4F_BINARY_DIR}/libsrc)
LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/fortran ${CMAKE_SOURCE_DIR}/libsrc
  ${NC4F_BINARY_DIR}/fortran)

# Create pre-compiled tests.
FOREACH(t ${TESTS})
  ADD_EXECUTABLE(${t} ${${t}_SOURCES})
  TARGET_LINK_LIBRARIES(${t} netcdff)
  ADD_TEST(${t} ${t})
ENDFOREACH()

# Copy test scripts in to current directory.
FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
FILE(COPY ${COPY_FILES}
	  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
	  FILE_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE)

# Add script-based tests.
FOREACH(t ${SCRIPT_TESTS})
  add_sh_test(nf_test ${t})
ENDFOREACH()


# test_get.F and test_put.f need to be distributed, so that the user
# need not have m4. fills.nc is used by test program ftest.
SET (EXTRA_DIST test_get.F test_put.F fills.cdl create_fills.sh
  run_f90_par_test.sh run_f77_par_test.sh ref_fills.nc)

# Cleaning up files created during the process.
SET (CLEANFILES scratch.nc test.nc copy.nc fills.nc tst_*.nc ftst_*.nc f90tst_*.nc)

SET (MAINTAINERCLEANFILES test_get.F test_put.F)
