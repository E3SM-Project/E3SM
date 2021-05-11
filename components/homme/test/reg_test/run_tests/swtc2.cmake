###############################################################
# RKSSP default benchmark (used to check nothing is broken)
###############################################################
#
# Spectral Element -- swtc2
# 
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME swtc2)
# The specifically compiled executable that this test uses
SET(EXEC_NAME swtcA)

SET(NUM_CPUS 16)

SET(MESH_FILES ${HOMME_ROOT}/test/mesh_files/mountain_10_x2.g)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)

SET(NC_OUTPUT_FILES swtc21.nc)
