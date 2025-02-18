# Build the CF validator tool

set (CF_STANDARD_NAME_FILE "${CMAKE_CURRENT_SOURCE_DIR}/cf-standard-name-table.yaml"
     CACHE STRING "Location of the cf standard name yaml table")
set (CF_SCREAM_NAME_FILE "${CMAKE_CURRENT_SOURCE_DIR}/cf-eamxx-name-table.yaml"
     CACHE STRING "Location of the scream-specific cf name yaml table")

add_executable(query-cf-database query-cf-database.cpp)
target_compile_definitions(query-cf-database PUBLIC
                           CF_STANDARD_NAME_FILE=${CF_STANDARD_NAME_FILE}
                           CF_SCREAM_NAME_FILE=${CF_SCREAM_NAME_FILE})

find_package (ekat HINTS ${EKAT_ROOT})
find_package (yaml-cpp HINTS ${YAML_CPP_ROOT})
target_link_libraries(query-cf-database ekat yaml-cpp)
