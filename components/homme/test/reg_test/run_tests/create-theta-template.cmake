# The name of this test (should be the basename of this file)
#example of name theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu
#or              theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu-kokkos
#adding BB to each var to avoid unwanted substitutions

SET(TEST_NAME BBNAME)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-nlev72BBIFKOKKOS)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/theta.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/acme-72*.ascii)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES
  jw_baroclinic1.nc
  jw_baroclinic2.nc)

# Specify test options, used to replace the cmake variables in the namelist
#DO NOT MOD
SET (HOMME_TEST_LIM 9)
SET (HOMME_TEST_MOISTURE dry)
SET (HOMME_TEST_HVSCALING 0) #const HV for now, tensor is tested in preqx

#mod
SET (HOMME_THETA_FORM BBADVFORM)
SET (HOMME_TTYPE BBTTYPE)
SET (HOMME_TEST_HVS BBHVS)
SET (HOMME_TEST_HVS_TOM BBTOM)
SET (HOMME_TEST_RSPLIT BBRSPLIT)
SET (HOMME_TEST_QSIZE BBQSIZE)
SET (HOMME_TEST_QSPLIT BBQSPLIT)
SET (HOMME_TEST_NUTOP BBNTOP)  #not BBNUTOP !
SET (HOMME_THETA_HY_MODE BBHYMODE)  

#DO NOT MOD
SET (HOMME_TEST_TIME_STEP 600)
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
