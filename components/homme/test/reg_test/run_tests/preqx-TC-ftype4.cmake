###############################################################
# preqx model, RK5 + PIO_INTERP 
###############################################################
#
# 1 day of the DCMIP2016 TC test, at NE=8
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME preqx-TC-ftype4)
# The specifically compiled executable that this test uses
SET(EXEC_NAME preqx-nlev30-interp)


SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/preqx-TC-ftype4.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*30*)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES   dcmip2016_test21.nc)



