set(P3_SRCS
  p3_f90.cpp
  p3_ic_cases.cpp
  micro_p3_iso_c.f90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/micro_p3.F90
)

set(P3_HEADERS
  p3_f90.hpp
  p3_ic_cases.hpp
)

# link_directories(${SCREAM_TPL_LIBRARY_DIRS} ${SCREAM_LIBRARY_DIRS})
add_library(p3 ${P3_SRCS})
target_include_directories(p3 PUBLIC ${SCREAM_INCLUDE_DIRS} ${SCREAM_TPL_INCLUDE_DIRS})
set_target_properties(p3 PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
# target_link_libraries(p3 scream_share ${SCREAM_TPL_LIBRARIES})

add_subdirectory(tests)
