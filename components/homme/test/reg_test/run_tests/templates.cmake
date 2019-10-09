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
SET(TEST_NAME templates)
# The specifically compiled executable that this test uses
SET(EXEC_NAME baroCam)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES 
  ${HOMME_ROOT}/test/reg_test/namelists/template1.nl
  ${HOMME_ROOT}/test/reg_test/namelists/template2.nl
)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*26*)
SET(MESH_FILES ${HOMME_ROOT}/test/mesh_files/mountain_10_x2.g)

SET(NC_OUTPUT_FILES 
  template1-held_suarez01.nc 
  template2-held_suarez01.nc 
)


