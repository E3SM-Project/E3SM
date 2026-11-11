# Homme system tests: shallow water equations (sweqx).
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

homme_add_system_test (swtc1
  EXEC          swtcA
  RANKS         16
  NAMELISTS     swtc1.nl
  OUTPUTS       swtc11.nc)

homme_add_system_test (swtc2
  EXEC          swtcA
  RANKS         16
  NAMELISTS     swtc2.nl
  MESH_FILES    ${HOMME_TEST_DATA_DIR}/mesh_files/mountain_10_x2.g
  OUTPUTS       swtc21.nc)

homme_add_system_test (swtc5
  EXEC          swtcB
  RANKS         16
  NAMELISTS     swtc5.nl
  NCL_FILES     swtc5ref.ncl
  OUTPUTS       swtc51.nc)

homme_add_system_test (swtc6
  EXEC          swtcA
  RANKS         16
  NAMELISTS     swtc6-short.nl swtc6-short-exodus.nl
  MESH_FILES    ${HOMME_TEST_DATA_DIR}/mesh_files/TEMPEST_NE2.g
  OUTPUTS       swtc61.nc
  REF_OUTPUTS   swtc61.nc
  CHECK_OUTPUTS exodus-swtc61.nc
  REF_TOL       5E-14)

homme_add_system_test (sweqx-dbvor
  EXEC          swtcC
  RANKS         4
  NAMELISTS     sweqx-dbvor.nl
  OUTPUTS       planar_dbl_vrtx1.nc)
