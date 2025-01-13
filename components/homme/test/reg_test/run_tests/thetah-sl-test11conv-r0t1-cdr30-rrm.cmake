# The name of this test (should be the basename of this file)
SET(TEST_NAME thetah-sl-test11conv-r0t1-cdr30-rrm)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-l-nlev30)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetah-sl-test11conv-r0t1-cdr30-rrm.nl)

SET(MESH_FILES ${HOMME_ROOT}/test/mesh_files/mountain_10_x2.g)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES 
  dcmip2012_test1_3a_conv1.nc)
