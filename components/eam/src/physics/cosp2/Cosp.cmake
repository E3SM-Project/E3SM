# In order for flag properties to persist, this file needs to be included,
# not added as a subdirectory

set(COSP_SOURCES
   eam/src/physics/cosp2/cosp_kinds.F90
   eam/src/physics/cosp2/external/src/cosp_constants.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_cloudsat_interface.F90
   eam/src/physics/cosp2/local/cosp_config.F90
   eam/src/physics/cosp2/local/cosp.F90
   eam/src/physics/cosp2/external/src/cosp_stats.F90
   eam/src/physics/cosp2/external/src/simulator/quickbeam/quickbeam.F90
   eam/src/physics/cosp2/external/src/simulator/parasol/parasol.F90
   eam/src/physics/cosp2/external/src/simulator/actsim/lidar_simulator.F90
   eam/src/physics/cosp2/external/src/simulator/icarus/icarus.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_atlid_interface.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_calipso_interface.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_grLidar532_interface.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_isccp_interface.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_misr_interface.F90
   eam/src/physics/cosp2/external/src/simulator/MISR_simulator/MISR_simulator.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_modis_interface.F90
   eam/src/physics/cosp2/local/modis_simulator.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_rttov_interfaceSTUB.F90
   eam/src/physics/cosp2/external/src/simulator/rttov/cosp_rttovSTUB.F90
   eam/src/physics/cosp2/external/src/simulator/cosp_parasol_interface.F90
   eam/src/physics/cosp2/subcol/scops.F90
   eam/src/physics/cosp2/subcol/prec_scops.F90
   eam/src/physics/cosp2/optics/cosp_utils.F90
   eam/src/physics/cosp2/optics/cosp_optics.F90
   eam/src/physics/cosp2/optics/quickbeam_optics.F90
   eam/src/physics/cosp2/subcol/mo_rng.F90
   eam/src/physics/cosp2/cosp_errorHandling.F90
   eam/src/physics/cosp2/optics/array_lib.F90
   eam/src/physics/cosp2/optics/math_lib.F90
   eam/src/physics/cosp2/optics/optics_lib.F90
   eam/src/physics/cosp2/optics/mrgrnk.F90
)

list(APPEND SOURCES ${COSP_SOURCES})
