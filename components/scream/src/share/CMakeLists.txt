set(SHARE_SRC
  scream_config.cpp
  scream_session.cpp
  atm_process/atmosphere_process_group.cpp
  atm_process/atmosphere_process_dag.cpp
  field/field_alloc_prop.cpp
  field/field_identifier.cpp
  field/field_header.cpp
  field/field_layout.cpp
  field/field_tracking.cpp
  grid/se_grid.cpp
  grid/simple_grid.cpp
  util/scream_test_session.cpp
  util/scream_time_stamp.cpp
)

add_library(scream_share ${SHARE_SRC})
target_include_directories(scream_share PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(scream_share SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS})
target_link_libraries(scream_share PUBLIC ekat)
set_target_properties(scream_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
# link_directories(${SCREAM_TPL_LIBRARY_DIRS})
# target_link_libraries(scream_share PUBLIC ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
