###############################################################
# Explicit leapfrog default benchmark
###############################################################
#
# Spectral Element -- swtc5
# NE=30, dt=90, nu=1.5e15, limiter=0, filter_freq=0, NP=4
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME swtc5-mr)
# The type of run (preqx,sweqx,swdgx,etc.)
SET(TEST_TYPE sweqx)
# The specifically compiled executable that this test uses
SET(EXEC_NAME swtcB)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)
SET(MESH_FILES ${HOMME_ROOT}/test/mesh_refine/grids/mountain_10_x8.g)

SET(NC_OUTPUT_FILES swtc51.nc)
