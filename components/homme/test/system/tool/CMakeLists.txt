# Homme system tests: homme_tool.
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

homme_add_system_test (hommetool
  EXEC          tool-nlev26
  RANKS         16
  NAMELISTS
    tool-mappingfiles.nl
    tool-template1.nl
    tool-template2.nl
    tool-topooutput.nl
    tool-toposmooth-gll.nl
    tool-zinterpolate.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*26*
  MESH_FILES
    ${HOMME_TEST_DATA_DIR}/mesh_files/mountain_10_x2.g
    ${HOMME_TEST_DATA_DIR}/mesh_files/9x16_scrip.nc
  OUTPUTS
    map_ne4np4_to_9x16_scrip_bilin.nc
    template1-ne4np4_tmp1.nc
    template2-ne0np4_tmp1.nc
    phis-baroclinic1.nc
    phis-smoothed1.nc
    phis-smoothed1.interp1.nc)
