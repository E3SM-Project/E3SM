# Kokkos version of pre-existing F90 test.
SET(TEST_NAME thetanh-moist-bubble-sl-pg2-kokkos)
SET(EXEC_NAME theta-l-nlev20-native-kokkos)
SET(NUM_CPUS 16)
SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetanh-moist-bubble-sl-pg2.nl)
SET(NC_OUTPUT_FILES planar_rising_bubble_pg21.nc)
