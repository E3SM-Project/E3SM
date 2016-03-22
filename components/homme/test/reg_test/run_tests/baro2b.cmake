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
SET(TEST_NAME baro2b)
# The specifically compiled executable that this test uses
SET(EXEC_NAME baroC)


SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)
SET(REFSOLN_FILES ${HOMME_ROOT}/test/reg_test/ref_sol/T340ref.nc)

# compare all of these files against baselines:
SET(NC_OUTPUT_FILES 
  baro2b-asp_baroclinic1.nc 
  baro2b-asp_baroclinic2.nc)

# replace these slow running thread tests with thread tests in baroCamMoist:
#SET(OMP_SUB_TESTS true)
#SET(OMP_NUM_THREADS 4)
#SET(OMP_NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}-omp.nl)


# compare openMP output vs single threaded output
#SET(NC_OUTPUT_REF   
#baro2b-asp_baroclinic1.nc  
#baro2b-asp_baroclinic2.nc  
#)
#SET(NC_OUTPUT_CHECKREF    
#baro2b-omp-asp_baroclinic1.nc 
#baro2b-omp-asp_baroclinic2.nc 
#)

