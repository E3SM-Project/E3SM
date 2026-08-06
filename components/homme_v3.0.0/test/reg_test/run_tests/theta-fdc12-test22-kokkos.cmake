# The name of this test (should be the basename of this file)
#example of name theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu
#or              theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu-kokkos
#adding BB to each var to avoid unwanted substitutions

SET(TEST_NAME theta-fdc12-test22-kokkos)
# The specifically compiled executable that this test uses
SET(EXEC_NAME theta-nlev128-kokkos)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/theta-fdc12-test22.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/sab*-128.ascii)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES
  dcmip2012_test2_21.nc)


#DO NOT MOD
SET (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
