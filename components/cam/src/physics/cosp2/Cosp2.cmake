# In order for flag properties to persist, this file needs to be included,
# not added as a subdirectory

set(COSP_SOURCES
  cam/src/physics/cosp2/llnl/cosp_radar.F90
  cam/src/physics/cosp2/cosp_types.F90
  cam/src/physics/cosp2/cosp_constants.F90
  cam/src/physics/cosp2/cosp_simulator.F90
  cam/src/physics/cosp2/cosp_utils.F90
  cam/src/physics/icarus-scops/scops.f90
  cam/src/physics/icarus-scops/icarus.f90
  cam/src/physics/cosp2/llnl/prec_scops.f90
  cam/src/physics/cosp2/cosp.F90
  cam/src/physics/cosp2/cosp_stats.F90
  cam/src/physics/cosp2/llnl/pf_to_mr.f90
  cam/src/physics/cosp2/cosp_lidar.F90
  cam/src/physics/cosp2/quickbeam/radar_simulator_types.f90
  cam/src/physics/cosp2/quickbeam/zeff.f90
  cam/src/physics/cosp2/quickbeam/array_lib.f90
  cam/src/physics/cosp2/quickbeam/atmos_lib.f90
  cam/src/physics/cosp2/quickbeam/dsd.f90
  cam/src/physics/cosp2/quickbeam/calc_Re.f90
  cam/src/physics/cosp2/quickbeam/format_input.f90
  cam/src/physics/cosp2/quickbeam/gases.f90
  cam/src/physics/cosp2/quickbeam/scale_LUTs_io.f90
  cam/src/physics/cosp2/quickbeam/radar_simulator_init.f90
  cam/src/physics/cosp2/quickbeam/math_lib.f90
  cam/src/physics/cosp2/quickbeam/mrgrnk.f90
  cam/src/physics/cosp2/quickbeam/optics_lib.f90
  cam/src/physics/cosp2/quickbeam/radar_simulator.f90
  cam/src/physics/cosp2/actsim/lidar_simulator.F90
  cam/src/physics/cosp2/llnl/llnl_stats.F90
  cam/src/physics/cosp2/actsim/lmd_ipsl_stats.F90
  cam/src/physics/cosp2/cosp_isccp_simulator.F90
  cam/src/physics/cosp2/cosp_misr_simulator.F90
  cam/src/physics/cosp2/MISR_simulator/MISR_simulator.f90
  cam/src/physics/cosp2/cosp_modis_simulator.F90
  cam/src/physics/cosp2/MODIS_simulator/modis_simulator.F90
)

set(COSP_FIXED
  cam/src/physics/icarus-scops/scops.f90
  cam/src/physics/icarus-scops/icarus.f90
  cam/src/physics/cosp2/llnl/prec_scops.f90
  cam/src/physics/cosp2/llnl/pf_to_mr.f90
  cam/src/physics/cosp2/MISR_simulator/MISR_simulator.f90
)

set(COSP_NOAUTO
  cam/src/physics/cosp2/quickbeam/mrgrnk.f90
)

foreach(ITEM IN LISTS COSP_SOURCES)
  list(FIND COSP_FIXED ${ITEM} ITEM_IS_FIXED)
  list(FIND COSP_NOAUTO ${ITEM} ITEM_IS_NO_AUTO)

  if (ITEM IN_LIST COSP_FIXED)
    e3sm_add_flags("${ITEM}" "${FIXEDFLAGS}")
  else()
    e3sm_add_flags("${ITEM}" "${FREEFLAGS}")
  endif()

  if (NOT ITEM IN_LIST COSP_NOAUTO)
    e3sm_add_flags("${ITEM}" "${FC_AUTO_R8}")
  endif()
endforeach()

list(APPEND SOURCES ${COSP_SOURCES})
