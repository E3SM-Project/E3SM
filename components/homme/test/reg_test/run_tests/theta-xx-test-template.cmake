# The name of this test (should be the basename of this file)
#theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu
SET(TEST_NAME BBNAME)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-nlev72)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/preqx.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/acme-72*.ascii)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES
  jw_baroclinic1.nc
  jw_baroclinic2.nc)

# Specify test options, used to replace the cmake variables in the namelist
#DO NOT MOD
SET (HOMME_TEST_LIM 9)
SET (HOMME_TEST_MOISTURE dry)

#mod
SET (HOMME_THETA_FORM ADVFORM)
SET (HOMME_TTYPE BBTTYPE)
SET (HOMME_TEST_HVSCALING AAHVS)
SET (HOMME_TEST_HVS_TOM AAHVST)
SET (HOMME_TEST_RSPLIT 3)
SET (HOMME_TEST_QSIZE 1)
SET (HOMME_TEST_NUTOP 0)
SET (HOMME_TEST_NU 1e15)
SET (HOMME_TEST_NUDIV 1e15)

#DO NOT MOD
SET (HOMME_TEST_TIME_STEP 600)
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
