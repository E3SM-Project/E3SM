###############################################################
# swim test case 5 (implicit regression test)
# np=4, ne=16, ndays=1, tstep=1200, tstep_type=13(BDF2)
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME swimtc5)
# The type of run (preqx,sweqx,swdgx,etc.)
SET(TEST_TYPE swim)
# The specifically compiled executable that this test uses
SET(EXEC_NAME swim5)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)
SET(TRILINOS_XML_FILE  ${HOMME_ROOT}/test/reg_test/trilinos_xml_files/swim/trilinosOptions.xml)
SET(NC_OUTPUT_FILES swtc51.nc)
