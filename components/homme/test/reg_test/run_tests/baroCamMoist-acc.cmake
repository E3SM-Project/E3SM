###############################################################
# RK + PIO_INTERP 
###############################################################
#
# Spectral Element -- 9 days of ASP baroclinic test
# (Jablonowski and Williamson test + 4 tracers)
# NE=15, dt=150, nu=1e16, filter_freq=0, NV=4, PLEV=26
# (explicit RK with subcycling)
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME baroCamMoist-acc)

# The specifically compiled executable that this test uses
SET(EXEC_NAME baroCam-acc)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES 
${HOMME_ROOT}/test/reg_test/namelists/baroCamMoist.nl
)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)

SET(NC_OUTPUT_FILES 
  camBaroMoist-asp_baroclinic1.nc 
  camBaroMoist-asp_baroclinic2.nc
)


# OMP tests (disabled for now)
#SET(OMP_SUB_TESTS true)
#SET(OMP_NUM_THREADS 4)
#SET(OMP_NAMELIST_FILES 
#${HOMME_ROOT}/test/reg_test/namelists/baroCamMoist-omp4.nl
#)


# compare output with CPU-only baroCamMoist test:
SET(TESTCASE_REF_TOL 1E-11)
SET(NC_OUTPUT_REF   
  ../../baroCamMoist/movies/camBaroMoist-asp_baroclinic1.nc 
  ../../baroCamMoist/movies/camBaroMoist-asp_baroclinic2.nc 
)
SET(NC_OUTPUT_CHECKREF    
  camBaroMoist-asp_baroclinic1.nc 
  camBaroMoist-asp_baroclinic2.nc
)


