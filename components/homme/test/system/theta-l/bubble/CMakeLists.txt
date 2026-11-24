# Homme system tests: theta-l rising bubble tests (planar).
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

# F90 and C++ versions. In BFB builds, also check that they give the same answers
homme_add_f90_cxx_system_tests (thetanh-moist-bubble CXX_VS_F90
  EXEC          theta-l-nlev20-native
  RANKS         4
  NAMELISTS     thetanh-moist-bubble.nl
  OUTPUTS       planar_rising_bubble1.nc)

homme_add_f90_cxx_system_tests (thetanh-dry-bubble CXX_VS_F90
  EXEC          theta-l-nlev20-native
  RANKS         4
  NAMELISTS     thetanh-dry-bubble.nl
  OUTPUTS       planar_rising_bubble1.nc)

# With semi-Lagrangian transport
if (HOMME_ENABLE_COMPOSE)
  homme_add_f90_cxx_system_tests (thetanh-moist-bubble-sl CXX_VS_F90
    EXEC          theta-l-nlev20-native
    RANKS         16
    NAMELISTS     thetanh-moist-bubble-sl.nl
    OUTPUTS       planar_rising_bubble1.nc)

  homme_add_f90_cxx_system_tests (thetanh-moist-bubble-sl-pg2 CXX_VS_F90
    EXEC          theta-l-nlev20-native
    RANKS         16
    NAMELISTS     thetanh-moist-bubble-sl-pg2.nl
    OUTPUTS       planar_rising_bubble_pg21.nc)
endif()
