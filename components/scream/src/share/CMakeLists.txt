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
  grid/point_grid.cpp
  util/scream_test_session.cpp
  util/scream_time_stamp.cpp
)

add_library(scream_share ${SHARE_SRC})
target_include_directories(scream_share PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(scream_share SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS})
target_link_libraries(scream_share PUBLIC ekat)
set_target_properties(scream_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})

# The "build_cf_dictionary" target downloads an XML file containing valid field
# names and aliases based on the CF conventions, and transforms it into a YAML
# file that we use to validate our own fields against it. This target is only
# available if CMake finds wget on your system.
find_program(WGET wget)
if (NOT WGET STREQUAL "WGET-NOTFOUND")
  set(CF_CONVENTIONS_XML_FILE cf-standard-name-table.xml)
  set(CF_CONVENTIONS_URL "http://cfconventions.org/Data/cf-standard-names/76/src/${CF_CONVENTIONS_XML_FILE}")
  add_custom_target(build_cf_dictionary
    ${WGET} ${CF_CONVENTIONS_URL}
    COMMAND ${PROJECT_SOURCE_DIR}/scripts/cf-xml-to-yaml ${CF_CONVENTIONS_XML_FILE}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data
    COMMENT "Downloading CF conventions XML file and building YAML dictionary...")
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
