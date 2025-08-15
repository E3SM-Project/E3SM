# Build cosp external library and eamxx_cosp interface for EAMxx

# Sources for EAMxx interface to COSP
set(COSP_SRCS
  eamxx_cosp.cpp
  cosp_c2f.F90
)
set(COSP_HEADERS
  eamxx_cosp.hpp
  cosp_functions.hpp
)

# Build external COSP library (this is all fortran code)
set(EAM_COSP_DIR ${SCREAM_BASE_DIR}/../eam/src/physics/cosp2)
set(EXTERNAL_SRC
   ${EAM_COSP_DIR}/cosp_kinds.F90
   ${EAM_COSP_DIR}/external/src/cosp_constants.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_cloudsat_interface.F90
   ${EAM_COSP_DIR}/external/src/cosp_config.F90
   ${EAM_COSP_DIR}/external/src/cosp.F90
   ${EAM_COSP_DIR}/external/src/cosp_stats.F90
   ${EAM_COSP_DIR}/external/src/simulator/quickbeam/quickbeam.F90
   ${EAM_COSP_DIR}/external/src/simulator/parasol/parasol.F90
   ${EAM_COSP_DIR}/external/src/simulator/actsim/lidar_simulator.F90
   ${EAM_COSP_DIR}/external/src/simulator/icarus/icarus.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_atlid_interface.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_calipso_interface.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_grLidar532_interface.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_isccp_interface.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_misr_interface.F90
   ${EAM_COSP_DIR}/external/src/simulator/MISR_simulator/MISR_simulator.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_modis_interface.F90
   ${EAM_COSP_DIR}/external/src/simulator/MODIS_simulator/modis_simulator.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_rttov_interfaceSTUB.F90
   ${EAM_COSP_DIR}/external/src/simulator/rttov/cosp_rttovSTUB.F90
   ${EAM_COSP_DIR}/external/src/simulator/cosp_parasol_interface.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/subcol/scops.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/subcol/prec_scops.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/optics/cosp_utils.F90 #optics/cosp_utils.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/optics/cosp_optics.F90 #optics/cosp_optics.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/optics/quickbeam_optics/quickbeam_optics.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/subcol/rng/mo_rng.F90
   ${EAM_COSP_DIR}/cosp_errorHandling.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/optics/quickbeam_optics/array_lib.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/optics/quickbeam_optics/math_lib.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/optics/quickbeam_optics/optics_lib.F90
   ${EAM_COSP_DIR}/external/subsample_and_optics_example/optics/quickbeam_optics/mrgrnk.F90
)
add_library(cosp ${EXTERNAL_SRC})

# Build interface code
add_library(eamxx_cosp ${COSP_SRCS})
target_link_libraries(eamxx_cosp physics_share scream_share cosp)
target_compile_options(eamxx_cosp PUBLIC)
target_compile_definitions(eamxx_cosp PUBLIC EAMXX_HAS_COSP)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE eamxx_cosp)
endif()
