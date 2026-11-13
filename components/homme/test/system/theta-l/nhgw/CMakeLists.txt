# Homme system tests: theta-l nonhydrostatic gravity waves (planar).
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

# The C++ versions are only tested in BFB builds, to compare against F90
homme_add_f90_cxx_system_tests (thetah-nhgw CXX_IF_BFB CXX_VS_F90
  EXEC          theta-l-nlev10-native
  RANKS         4
  NAMELISTS     thetah-nhgw.nl
  OUTPUTS       planar_nonhydro_gravity_wave1.nc)

homme_add_f90_cxx_system_tests (thetanh-nhgw CXX_IF_BFB CXX_VS_F90
  EXEC          theta-l-nlev10-native
  RANKS         4
  NAMELISTS     thetanh-nhgw.nl
  OUTPUTS       planar_nonhydro_gravity_wave1.nc)

homme_add_f90_cxx_system_tests (thetah-nhgw-slice CXX_IF_BFB CXX_VS_F90
  EXEC          theta-l-nlev10-native
  RANKS         4
  NAMELISTS     thetah-nhgw-slice.nl
  OUTPUTS       planar_nonhydro_gravity_wave1.nc)

homme_add_f90_cxx_system_tests (thetanh-nhgw-slice CXX_IF_BFB CXX_VS_F90
  EXEC          theta-l-nlev10-native
  RANKS         4
  NAMELISTS     thetanh-nhgw-slice.nl
  OUTPUTS       planar_nonhydro_gravity_wave1.nc)
