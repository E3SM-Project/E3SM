INCLUDE (ScreamUtils)

# Create the test
SET (NEED_LIBS shoc cld_fraction p3 scream_rrtmgp scream_control scream_share physics_share)
CreateUnitTest(shoc_cld_p3_rrtmgp shoc_cld_p3_rrtmgp.cpp "${NEED_LIBS}" LABELS "shoc;cld;p3;rrtmgp;physics")

# Set AD configurable options
SetVarDependingOnTestProfile(NUM_STEPS 2 24 48)
set (ATM_TIME_STEP 1800)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)

## Copy input files to run directory
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
set (RRTMGP_EXT_DIR ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shoc_cld_p3_rrtmgp_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_cld_p3_rrtmgp_init_ne2np4.nc COPYONLY)

configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)

configure_file(${RRTMGP_EXT_DIR}/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/examples/all-sky/rrtmgp-allsky.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
