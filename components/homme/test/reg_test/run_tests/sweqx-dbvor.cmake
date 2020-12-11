###############################################################
#
# 250 steps of the low-resolution (20x20 elements) planar double vortex test case
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME sweqx-dbvor)
# The specifically compiled executable that this test uses
SET(EXEC_NAME swtcC)

SET(NUM_CPUS 4)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/sweqx-dbvor.nl)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES planar_dbl_vrtx1.nc)


