###############################################################
# Explicit leapfrog default benchmark
###############################################################
#
# swtc6 (Rossby Wave)
# short time, low res, check conservation

# The name of this test (should be the basename of this file)
SET(TEST_NAME swtc6)
# The type of run (preqx,sweqx,swdgx,etc.)
SET(TEST_TYPE sweqx)
# The specifically compiled executable that this test uses
SET(EXEC_NAME swtcC)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/swtc6-short.nl)

SET(NC_OUTPUT_FILES swtc61.nc)

