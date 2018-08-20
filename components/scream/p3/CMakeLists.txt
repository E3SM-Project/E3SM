set(P3_SRCS
  p3_f90.cpp
  p3_ic_cases.cpp
  micro_p3_iso_c.f90
  ../../cam/src/physics/cam/micro_p3.F90
)

set(P3_HEADERS
  p3_f90.hpp
  p3_ic_cases.hpp
)

add_library(p3 ${P3_SRCS})
target_include_directories(p3 PUBLIC ${SCREAM_INCLUDES})
set_target_properties(p3 PROPERTIES
  LINK_FLAGS ${SCREAM_LINK_FLAGS}
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(p3 scream_util ${SCREAM_LIBRARIES})

add_subdirectory(tests)
