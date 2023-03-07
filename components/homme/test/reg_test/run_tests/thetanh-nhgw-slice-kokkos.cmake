# Kokkos version of pre-existing F90 test.
SET(TEST_NAME thetanh-nhgw-slice-kokkos)
SET(EXEC_NAME theta-l-nlev10-native-kokkos)
SET(NUM_CPUS 4)
SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetanh-nhgw-slice.nl)
SET(NC_OUTPUT_FILES planar_nonhydro_gravity_wave1.nc)
