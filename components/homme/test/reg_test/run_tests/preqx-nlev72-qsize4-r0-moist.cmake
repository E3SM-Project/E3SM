# The name of this test (should be the basename of this file)
SET(TEST_NAME preqx-nlev72-qsize4-r0-moist)
# The specifically compiled executable that this test uses
SET(EXEC_NAME preqx-nlev72)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/preqx.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/acme-72*)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES
  jw_baroclinic1.nc
  jw_baroclinic2.nc)

# Specify test options, used to replace the cmake variables in the namelist
SET (HOMME_TEST_QSIZE 4)
SET (HOMME_TEST_RSPLIT 0)
SET (HOMME_TEST_MOISTURE notdry)
SET (HOMME_TEST_TIME_STEP 300)
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
