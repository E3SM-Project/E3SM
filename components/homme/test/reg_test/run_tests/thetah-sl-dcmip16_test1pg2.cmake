# The name of this test (should be the basename of this file)
SET(TEST_NAME thetah-sl-dcmip16_test1pg2)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-l-nlev30)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetah-sl-dcmip16_test1pg2.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*30*)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES 
  dcmip2016_test1_pg21.nc)
