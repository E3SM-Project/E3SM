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
# The type of run (preqx,sweqx,swdgx,etc.)
SET(TEST_TYPE preqx)
# The specifically compiled executable that this test uses
SET(EXEC_NAME baroC)


SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/${NAMELIST_DIR}/${TEST_NAME}.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)
SET(REFSOLN_FILES ${HOMME_ROOT}/test/reg_test/ref_sol/T340ref.nc)

SET(OMP_SUB_TESTS true)
SET(OMP_NUM_THREADS 4)
SET(OMP_NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/${NAMELIST_DIR}/${TEST_NAME}-omp.nl)

SET(NC_OUTPUT_FILES asp_baroclinic1.nc asp_baroclinic2.nc)
