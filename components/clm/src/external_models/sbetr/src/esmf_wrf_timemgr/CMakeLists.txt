set(BETR_ESMF_SOURCES
  ESMF.F90
  ESMF_AlarmClockMod.F90
  ESMF_AlarmMod.F90
  ESMF_BaseMod.F90
  ESMF_BaseTimeMod.F90
  ESMF_CalendarMod.F90
  ESMF_ClockMod.F90
  ESMF_FractionMod.F90
  ESMF_ShrTimeMod.F90
  ESMF_Stubs.F90
  ESMF_TimeIntervalMod.F90
  ESMF_TimeMod.F90
  MeatMod.F90
  wrf_error_fatal.F90
  wrf_message.F90
  )

configure_file(ESMF_TimeMgr.inc ESMF_TimeMgr.inc COPYONLY)
configure_file(ESMF_Macros.inc ESMF_Macros.inc COPYONLY)
include_directories(${CMAKE_BINARY_DIR}/src/esmf_wrf_timemgr)
include_directories(${CMAKE_BINARY_DIR}/src/shr)

include(add_betr_library)
add_betr_library(betr_esmf ${BETR_ESMF_SOURCES})

# disable all compiler warnings for clm code
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

set(BETR_LIBRARIES betr_esmf;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_esmf;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_esmf DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/esmf)
endif()
