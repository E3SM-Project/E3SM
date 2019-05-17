set(CMAKE_Fortran_FLAGS "${FFLAGS} ${INCLDIR} -I${CMAKE_BINARY_DIR}")
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_Fortran_COMPILER ${FC})

set(COSP_SOURCES
  llnl/cosp_radar.F90
  cosp_types.F90
  cosp_constants.F90
  cosp_simulator.F90
  cosp_utils.F90
  ../icarus-scops/scops.f90
  ../icarus-scops/icarus.f90
  llnl/prec_scops.f90
  cosp.F90
  cosp_stats.F90
  llnl/pf_to_mr.f90
  cosp_lidar.F90
  quickbeam/radar_simulator_types.f90
  quickbeam/zeff.f90
  quickbeam/array_lib.f90
  quickbeam/atmos_lib.f90
  quickbeam/dsd.f90
  quickbeam/calc_Re.f90
  quickbeam/format_input.f90
  quickbeam/gases.f90
  quickbeam/scale_LUTs_io.f90
  quickbeam/radar_simulator_init.f90
  quickbeam/math_lib.f90
  quickbeam/mrgrnk.f90
  quickbeam/optics_lib.f90
  quickbeam/radar_simulator.f90
  actsim/lidar_simulator.F90
  llnl/llnl_stats.F90
  actsim/lmd_ipsl_stats.F90
  cosp_isccp_simulator.F90
  cosp_misr_simulator.F90
  MISR_simulator/MISR_simulator.f90
  cosp_modis_simulator.F90
  MODIS_simulator/modis_simulator.F90
)

set(COSP_FIXED
  ../icarus-scops/scops.f90
  ../icarus-scops/icarus.f90
  llnl/prec_scops.f90
  llnl/pf_to_mr.f90
  MISR_simulator/MISR_simulator.f90
)

set(COSP_NOAUTO
  quickbeam/mrgrnk.f90
)

foreach(ITEM IN LISTS COSP_SOURCES)
  list(FIND COSP_FIXED ${ITEM} ITEM_IS_FIXED)
  list(FIND COSP_NOAUTO ${ITEM} ITEM_IS_NO_AUTO)

  if (ITEM IN_LIST COSP_FIXED)
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FIXEDFLAGS}")
  else()
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FREEFLAGS}")
  endif()

  if (NOT ITEM IN_LIST COSP_NOAUTO)
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FC_AUTO_R8}")
  endif()
endforeach()

add_library(cosp)
target_sources(cosp PRIVATE ${COSP_SOURCES})
separate_arguments(FPPDEFS_LIST UNIX_COMMAND "${FPPDEFS}")
target_compile_definitions(cosp PRIVATE ${FPPDEFS_LIST})
