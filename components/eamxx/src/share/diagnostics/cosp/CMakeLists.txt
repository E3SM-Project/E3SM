# The COSP satellite simulator, as a diagnostic.
#
# COSP produces a whole suite of fields from one (expensive) pass, so CospDiag is
# a multi-output diagnostic: a request for any of the fields it computes builds
# this one diag, and they are all filled together. See eamxx_cosp_diag.hpp.
#
# This is a library of its own rather than part of eamxx_diagnostics, so that
# eamxx_diagnostics -- which everything under share/ links -- does not acquire a
# dependency on COSP's fortran. share/CMakeLists.txt hangs it off scream_share,
# which is what makes EAMXX_HAS_COSP_DIAG visible to register_diagnostics().
#
# NOTE: this builds at any precision. COSP's fortran interface is unconditionally
#       double (see cosp_functions.hpp), and the C++ side stages through double
#       buffers to match, so there is no SCREAM_DOUBLE_PRECISION gate here.

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

# The diagnostic and its C++/fortran bridge
add_library(eamxx_cosp_diag
  eamxx_cosp_diag.cpp
  cosp_c2f.F90
)

target_link_libraries(eamxx_cosp_diag PUBLIC
  eamxx_diagnostics
  eamxx_physics_share
  cosp)

target_compile_definitions(eamxx_cosp_diag PUBLIC EAMXX_HAS_COSP_DIAG)
