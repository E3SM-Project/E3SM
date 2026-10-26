set(TEST_NAME planar-transport-a-etm-kokkos)
set(EXEC_NAME theta-l-nlev128-native-kokkos)
set(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/planar-transport-a-etm.nl)
set(NUM_CPUS 16)
# 640 elements (ne_x=128, ne_y=5) with semi_lagrange_halo=3 requires >= 4 elements/rank
set(MAX_NUM_CPUS 160)

set(NC_OUTPUT_FILES planar_transport_a1.nc)

check_transport_error_norms(
  ${TEST_NAME} "planar_conv" "${TEST_NAME}_1.out"
  "3.9e-3;3.9e-3;4.1e-3")
