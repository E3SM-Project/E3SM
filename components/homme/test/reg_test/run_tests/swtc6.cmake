###############################################################
# Explicit leapfrog default benchmark
###############################################################
#
# swtc6 (Rossby Wave)
# short time, low res, check conservation

# The name of this test (should be the basename of this file)
SET(TEST_NAME swtc6)
# The specifically compiled executable that this test uses
SET(EXEC_NAME swtcA)

SET(NUM_CPUS 16)

SET(MESH_FILES ${HOMME_ROOT}/test/mesh_files/TEMPEST_NE2.g)

SET(NAMELIST_FILES 
   ${HOMME_ROOT}/test/reg_test/namelists/swtc6-short.nl
   ${HOMME_ROOT}/test/reg_test/namelists/swtc6-short-exodus.nl
)

# files to compare against baselines
SET(NC_OUTPUT_FILES swtc61.nc )

# compare exodus output vs internal cubed sphere output:
SET(TESTCASE_REF_TOL 1E-14)
#SET(TESTCASE_REF_TOL 1E-17)
SET(NC_OUTPUT_REF swtc61.nc )
SET(NC_OUTPUT_CHECKREF  exodus-swtc61.nc)


