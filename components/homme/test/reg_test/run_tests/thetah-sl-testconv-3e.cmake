# The name of this test (should be the basename of this file)
SET(TEST_NAME thetah-sl-testconv-3e)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-l-nlev30)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetah-sl-testconv-3e.nl)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES dcmip2012_test1_3e_conv1.nc)

check_transport_error_norms(
  ${TEST_NAME} "test1_conv" "${TEST_NAME}_1.out"
  "9.1e-2;1.1e-1;4.0e-2;8.1e-2")
 