# The name of this test (should be the basename of this file)
SET(TEST_NAME thetah-sl-test11conv-r1t2-cdr20)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-l-nlev30)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetah-sl-test11conv-r1t2-cdr20.nl)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES 
  dcmip2012_test1_3a_conv1.nc)
