
# The name of this test (should be the basename of this file)
SET(TEST_NAME thetaplanar-nlev128)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-l-nlev128-native)

SET(NUM_CPUS 4)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/theta-p3.nl)
SET(P3_FILES ${HOMME_ROOT}/p3tables/*)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES planar1.nc)


