set(TEST_NAME planar-transport-a-etm)
set(EXEC_NAME theta-l-nlev128-native)
set(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/planar-transport-a-etm.nl)
set(NUM_CPUS 16)

set(NC_OUTPUT_FILES planar_transport_a1.nc)

check_transport_error_norms(
  ${TEST_NAME} "planar_conv" "${TEST_NAME}_1.out"
  "3.9e-3;3.9e-3;4.1e-3")
