# Kokkos version of pre-existing F90 test.
SET(TEST_NAME thetanh-moist-bubble-kokkos)
SET(EXEC_NAME theta-l-nlev20-native-kokkos)
SET(NUM_CPUS 4)
SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/thetanh-moist-bubble.nl)
SET(NC_OUTPUT_FILES planar_rising_bubble1.nc)
