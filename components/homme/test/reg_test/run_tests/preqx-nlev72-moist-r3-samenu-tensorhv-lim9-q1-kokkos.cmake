# The name of this test (should be the basename of this file)
SET(TEST_NAME preqx-nlev72-moist-r3-samenu-tensorhv-lim9-q1-kokkos)
# The specifically compiled executable that this test uses
SET(EXEC_NAME preqx-nlev72-kokkos)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/preqx.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/acme-72*.ascii)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES
  jw_baroclinic1.nc
  jw_baroclinic2.nc)

# Specify test options, used to replace the cmake variables in the namelist
SET (HOMME_TEST_LIM 9)
SET (HOMME_TEST_QSIZE 1)
SET (HOMME_TEST_RSPLIT 3)
SET (HOMME_TEST_MOISTURE moist)

#const HV is HVSCALING=0, tensor HV is 3.2
#const HV: all preqxx tests use nu=7e15, to test nu!=nudiv, set nudiv to 1e15
#tensor HV: nu=nudiv=1e-9
SET (HOMME_TEST_HVSCALING 3.2)
SET (HOMME_TEST_NU 4.9e-7)
SET (HOMME_TEST_NUDIV 4.9e-7)
SET (HOMME_TEST_NUTOP 0)

SET (HOMME_TEST_TIME_STEP 300)
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
