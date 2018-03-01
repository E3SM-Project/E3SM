set(BETR_MATH_SOURCES
  func_data_type_mod.F90
  FindRootMod.F90
  InterpolationMod.F90
  MathfuncMod.F90
  ODEMod.F90
  BTridiagonalMod.F90
  KineticsMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)


include(add_betr_library)
add_betr_library(betr_math ${BETR_MATH_SOURCES})

set(BETR_LIBRARIES betr_math;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_math;${BETR_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_math DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/math)
endif()
