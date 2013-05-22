###############################################################
# Restart, and energy conservation test
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME baro2c)
# The type of run (preqx,sweqx,swdgx,etc.)
SET(TEST_TYPE preqx)
# The specifically compiled executable that this test uses
SET(EXEC_NAME baroC)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES
  ${HOMME_ROOT}/test/reg_test/${NAMELIST_DIR}/${TEST_NAME}-run1.nl
  ${HOMME_ROOT}/test/reg_test/${NAMELIST_DIR}/${TEST_NAME}-run2.nl
)

SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)

SET(OMP_SUB_TESTS true)
SET(OMP_NUM_THREADS 4)
SET(OMP_NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/${NAMELIST_DIR}/${TEST_NAME}-run2-omp.nl)

