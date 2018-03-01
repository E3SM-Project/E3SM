set(BETR_UTIL_SOURCES
  gbetrType.F90
  betr_constants.F90
  babortutils.F90
  BetrStatusType.F90
  BeTR_decompMod.F90
  betr_varcon.F90
  betr_ctrl.F90
  bshr_const_mod.F90
  bshr_infnan_mod.F90
  bshr_kind_mod.F90
  bshr_log_mod.F90
  bshr_nl_mod.F90
  bshr_string_mod.F90
  bshr_mpi_mod.F90
  bshr_assert_mod.F90
  bshr_sys_mod.F90
  bshr_timer_mod.F90
  betr_time_mod.F90
  betr_regression_mod.F90
  betr_utils.F90
  Tracer_varcon.F90
  )

configure_file(bshr_assert.h bshr_assert.h COPYONLY)
configure_file(bmpif.h bmpif.h COPYONLY)
configure_file(dtypes.h dtypes.h COPYONLY)

include(add_betr_library)
add_betr_library(betr_util ${BETR_UTIL_SOURCES})

set(BETR_LIBRARIES betr_util;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_util;${BETR_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/util)
endif()
