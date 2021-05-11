###############################################################
# RK + PIO (native grid output)
###############################################################
#
# Test the HOMME grid template generation code
# This runs 1 timstep, and then outputs details about the grid
# and the GLL dual grid 
# 
# Not tested:  From this output, NCL scripts are used to convert
# the data into a "*_latlon.nc" file used by CAM's interpic_new utility,
# and a "*_scip.nc" file used by ESMF and Tempest utilities that
# generate mapping files.
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME hommetool)
# The specifically compiled executable that this test uses
SET(EXEC_NAME tool-nlev26)

SET(NUM_CPUS 16)

# run in alphebetical order:
SET(NAMELIST_FILES 
  ${HOMME_ROOT}/test/reg_test/namelists/tool-mappingfiles.nl
  ${HOMME_ROOT}/test/reg_test/namelists/tool-template1.nl
  ${HOMME_ROOT}/test/reg_test/namelists/tool-template2.nl
  ${HOMME_ROOT}/test/reg_test/namelists/tool-topooutput.nl
  ${HOMME_ROOT}/test/reg_test/namelists/tool-toposmooth-gll.nl
  ${HOMME_ROOT}/test/reg_test/namelists/tool-zinterpolate.nl
)

SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)
SET(MESH_FILES ${HOMME_ROOT}/test/mesh_files/mountain_10_x2.g ${HOMME_ROOT}/test/mesh_files/9x16_scrip.nc)

SET(NC_OUTPUT_FILES 
  map_ne4np4_to_9x16_scrip_bilin.nc
  template1-ne4np4_tmp1.nc
  template2-ne0np4_tmp1.nc
  phis-baroclinic1.nc
  phis-smoothed1.nc
  phis-smoothed1.interp1.nc
)
