###############################################################
#
# Spectral Element -- 9 days of ASP baroclinic test
# (Jablonowski and Williamson test + 4 tracers)
# NE=15, dt=180, nu=1e16, filter_freq=0, NP=4, PLEV=26
# (explicit leap-frog)
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME baro2d-imp)

# The specifically compiled executable that this test uses
SET(EXEC_NAME primA)

# Number of CPUs for the test
SET(NUM_CPUS 16)

# The namelist file for the test
SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)

# The trilinos options file for the test
SET(TRILINOS_XML_FILE  ${HOMME_ROOT}/test/reg_test/trilinos_xml_files/prim/trilinosOptions.xml)

# The vertical coordiante files for the test
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)

# The output files produced from this run
SET(NC_OUTPUT_FILES baro2d-imp-asp_baroclinic1.nc baro2d-imp-asp_baroclinic2.nc)
