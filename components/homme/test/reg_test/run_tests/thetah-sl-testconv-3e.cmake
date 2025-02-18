# The name of this test (should be the basename of this file)
SET(TEST_NAME thetah-sl-testconv-3e)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-l-nlev30)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetah-sl-testconv-3e.nl)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES dcmip2012_test1_3e_conv1.nc)
