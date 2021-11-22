
SET(TEST_NAME theta-fhs1)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-nlev128)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetahs1.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/sab*-128.ascii)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES
  held_suarez01.nc
  held_suarez02.nc)

#DO NOT MOD
SET (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
