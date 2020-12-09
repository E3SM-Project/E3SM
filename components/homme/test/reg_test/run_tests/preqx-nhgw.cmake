###############################################################
#
# 200 steps of the low-resolution (33x33 elements) planar nonhydrostatic gravity wave test case
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME thetanh-nhgw)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-l-nlev10-native)

SET(NUM_CPUS 4)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetanh-nhgw.nl)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES planar_nonhydro_gravity_wave1.nc)


