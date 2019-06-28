set(EMI_SHR_SOURCES
  shr_assert_mod.F90
  shr_const_mod.F90
  shr_infnan_mod.F90
  shr_kind_mod.F90
  shr_log_mod.F90
  shr_sys_mod.F90
 )

include(add_emi_library)
add_emi_library(emi_shr ${EMI_SHR_SOURCES})

set(EMI_LIBRARIES emi_shr;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_shr;${EMI_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_shr DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include)
endif()

