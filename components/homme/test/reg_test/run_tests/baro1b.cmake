###############################################################
# Semi-implicit + PIO
###############################################################
#
# Spectral Element -- Polvani et al baroclinic test
# NE=9, dt=600, nu=7e5, filter_freq=0, NP=8, PLEV=20
# (semi-implicit)
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME baro1b)
# The specifically compiled executable that this test uses
SET(EXEC_NAME baroB)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*20*)
SET(NCL_FILES ${HOMME_ROOT}/test/reg_test/ncl/${TEST_NAME}.ncl)

SET(NC_OUTPUT_FILES baro1b-baroclinic1.nc)
