# The name of this test (should be the basename of this file)
#example of name theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu
#or              theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu-kokkos
#adding BB to each var to avoid unwanted substitutions

SET(TEST_NAME theta-f1-tt10-hvs1-hvst0-r3-qz1-nutopoff-kokkos)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-nlev128-kokkos)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/theta.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/sab*-128.ascii)

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
SET (HOMME_THETA_FORM 1)
SET (HOMME_TTYPE 10)
SET (HOMME_TEST_HVS 1)
SET (HOMME_TEST_HVS_TOM 0)
SET (HOMME_TEST_RSPLIT 3)
SET (HOMME_TEST_QSIZE 1)
SET (HOMME_TEST_QSPLIT )
SET (HOMME_TEST_NUTOP 0)  #not BBNUTOP !
SET (HOMME_THETA_HY_MODE false)  

#DO NOT MOD
SET (HOMME_TEST_TIME_STEP 600)
SET (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
