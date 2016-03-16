###############################################################
# Restart, and energy conservation test
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME baro2c)
# The specifically compiled executable that this test uses
SET(EXEC_NAME baroC)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES
  ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}-run1.nl
  ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}-run2.nl
)

SET(NC_OUTPUT_FILES 
  baro2c-run1-jw_baroclinic1.nc
  baro2c-run2-jw_baroclinic-000000000-1.nc)

SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)

SET(OMP_SUB_TESTS true)
SET(OMP_NUM_THREADS 4)
SET(OMP_NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}-run2-omp.nl)

#SET(OMP_NC_OUTPUT_FILES 
#  baro2c-run2-omp-jw_baroclinic-000000000-1.nc)


# compare openMP output vs single threaded output
SET(NC_OUTPUT_REF  baro2c-run2-jw_baroclinic-000000000-1.nc )
SET(NC_OUTPUT_CHECKREF baro2c-run2-omp-jw_baroclinic-000000000-1.nc )


